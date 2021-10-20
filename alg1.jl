using LinearAlgebra
using Laplacians
using SparseArrays
using Arpack
include("edgecore.jl")

function greedy(G,L,k,rep)
    x=argmax(L)[1];
    s=zeros(k);
    s[1]=x;
    n=G.n;
    sel=zeros(n);
    sel[x]=1;
    D=spzeros(n,n);
    D[x,x]=1;
    for i=2:k
        sco_V=zeros(n);
        e=eigcal(L+D);
		#e=fiedler(L+D)
		for j=1:n
			if sel[j]==1
				e[2][j]=0;
			end
		end
        for j=1:G.m
            u=G.u[j];v=G.v[j];
            sco_V[u]+=e[2][u]*e[2][v];
            sco_V[v]+=e[2][u]*e[2][v];
        end
        for j=1:G.n
            #sco_V[j]+=e[2][j]*e[2][j];
        end
        x=s[1];
        while sel[Int(x)]==1
            sco_V[Int(x)]=0;
            x=argmax(sco_V);
        end
        s[i]=x;
        sel[x]=1;
        D[x,x]=1;
    end
	if k>=2
		for i=1:rep
			s=update(G,L,k,s);
		end
	end
    return s;
end


function greedy2(G,L,k)
    x=argmax(L)[1];
    s=zeros(k);
    s[1]=x;
    n=G.n;
    sel=zeros(n);
    sel[x]=1;
    D=spzeros(n,n);
    D[x,x]=1;
    for i=2:k
        sco_V=zeros(n);
        e=eigcal(L+D);
        for j=1:G.n
            sco_V[j]+=e[2][j]^2;
        end
        x=s[1];
        while sel[Int(x)]==1
            sco_V[Int(x)]=0;
            x=argmax(sco_V);
        end
        s[i]=x;
        sel[x]=1;
        D[x,x]=1;
    end
    return s;
end

function update(G,L,k,s)
	n=G.n;
    snew=zeros(k);
    for i=1:k
        tmp=union(snew[1:i-1],s[i+1:k]);
        D=spzeros(n,n);
		sel=zeros(n);
        for j=1:i-1
			x=Int(snew[j]);
            D[x,x]=1;
			sel[x]=1;
        end
		for j=i+1:k
			x=Int(s[j]);
			D[x,x]=1;
			sel[x]=1;
		end
        sco_V=zeros(n);
        e=eigcal(L+D);
        for j=1:G.m
            u=G.u[j];v=G.v[j];
            sco_V[u]+=e[2][u]*e[2][v];
            sco_V[v]+=e[2][u]*e[2][v];
        end
		for j=1:n
			if sel[j]==1
				sco_V[j]=0;
			end
		end

        x=argmax(sco_V);
		D[x,x]=1;
		la=eigcal(L+D)[1][1];

		D[x,x]=0;
		y=Int(s[i]);
		D[y,y]=1;
		lo=eigcal(L+D)[1][1];
		if la>lo
			snew[i]=x;
		else
			snew[i]=y;
		end
    end
    return snew;
end

function dis(G,p)
    n=G.n;
    d=n*ones(n);
    d[p]=0;
    vis=zeros(n);
    vis[p]=1;
    v=1;
    b=zeros(n+1);
    b[1]=p;
    f=1;r=0;
    while f>r
        r+=1;
        x=Int(b[r]);
        for i=1:length(G.nbr[x])
            if vis[G.nbr[x][i]]==0
                f+=1;
                b[f]=G.nbr[x][i];
                d[G.nbr[x][i]]=d[x]+1;
                vis[G.nbr[x][i]]=1;
            end
        end
    end
    return d;
end

function kcent(G,L,k)
    n=G.n;
    d=n*ones(n);
    dnew=zeros(n);
    sel=zeros(k);
    x=argmax(L)[1];
    sel[1]=x;
    for i=2:k
        dnew=dis(G,Int(sel[i-1]))
        for j=1:n
            d[j]=min(d[j],dnew[j]);
            #d[j]+=dnew[j];
        end
        x=argmax(d);
        sel[i]=x;
    end
    return sel;
end


function ransam(p,n)
	s=sum(p);
	c=s*rand();
	tot=0;
	cho=n;
	for i=1:n
		tot+=p[i];
		if tot>=c
			cho=i;
			break;
		end
	end
	return cho;
end


function samp(G,L,k)
    n=G.n;
    d=zeros(n);
    s=zeros(k);
    for i=1:n
        d[i]=L[i,i];
    end
    #s[1]=ransam(d,n);
	s[1]=argmax(d);
    sel=zeros(n);
    sel[Int(s[1])]=1;
    S=union(1:n);
    setdiff!(S,s[1]);
    for i=2:k
		D=spzeros(n,n);
		for j=1:i-1
			x=Int(s[j]);
			D[x,x]=1;
		end
        Lg=L+D;
        e=eigcal(Lg);
        u=e[2];
        sco_V=zeros(n);


        for j=1:G.m
            x=G.u[j];y=G.v[j];
            sco_V[x]+=u[x]*u[y];
            sco_V[y]+=u[x]*u[y];
        end

        for j=1:n
            if sel[j]==1
                sco_V[j]=0;
            end
        end

        #x=Int(ransam(sco_V,n));
		x=argmax(sco_V);

        s[i]=x;
        sel[x]=1;
    end
	#println(s);
    re=k;
	succ=0;
	fail=0;
	upda=1;
	sco_V=zeros(n);
	del_V=zeros(n);
	e=zeros(n);
	D=spzeros(n,n);
    #for i=1:re
	while fail<100
		if upda==1
			#D=spzeros(n,n);
			for j=1:k
				x=Int(s[j]);
				D[x,x]=1;
			end
        	Lg=L+D;
        	e=eigcal(Lg);
        	u=e[2];
        	tot=0;

        	for j=1:G.m
            	sco_V[G.u[j]]+=u[G.u[j]]*u[G.v[j]];
            	sco_V[G.v[j]]+=u[G.u[j]]*u[G.v[j]];
        	end
        	for j=1:n
            	if sel[j]==1
					del_V[j]=1/sco_V[j];
                	sco_V[j]=0;
            	end
        	end
		end
        x=Int(ransam(sco_V,n));
        y=Int(ransam(del_V,n));
        lb=e[1][1];
		D[x,x]=1;D[y,y]=0;
        Lg=L+D;
        ee=eigcal(Lg)
        if lb<eigcal(Lg)[1][1]
            sel[x]=1;
            sel[y]=0;
            for j=1:k
                if Int(s[j])==y
                    s[j]=x;
                    #break;
                end
            end
			succ+=1;
			upda=1;
        else
			D[x,x]=0;D[y,y]=1;
			fail+=1;
			upda=0;
        end

    end
	#fout = open("ans.txt", "a")
	#println(fout);
	println(k,' ',succ,' ',fail)
	#close(fout);
    return s;
end


function fastgreedy(G,L,k,eps)
	x=argmax(L)[1];
	s=zeros(k);
	n=G.n;
	sel=zeros(n);
	sel[x]=1;
	s[1]=x;
	D=spzeros(n,n);
	D[x,x]=1;
	u=eigcal(L+D)[2];
	sco_V=zeros(n);
	for j=1:G.m
		sco_V[G.u[j]]+=u[G.u[j]]*u[G.v[j]];
		sco_V[G.v[j]]+=u[G.u[j]]*u[G.v[j]];
	end
	sco_V[x]=0;
	xx=argmax(sco_V);
	d=sco_V[xx];
	tot=2;
	s[2]=xx;
	sel[xx]=1;
	sco_V[xx]=0;
	while tot<k
		maxx=0;
		ss=0;
		rad=Int(round(rand()*(G.n-1)+1))
		for iii=1:n
			i=iii+rad;
			if i>n
				i=i-n;
			end
			if sco_V[i]>maxx
				maxx=sco_V[i];
			end
			if sco_V[i]>=d && tot<k
				tot+=1;
				sel[i]=1;
				s[tot]=i;
				ss=1;
			end
		end
		if ss==1
			d=d*(1-eps);
			D=spzeros(n,n);
			for i=1:tot
				x=Int(s[i]);
				D[x,x]+=1;
			end
			u=eigcal(D+L)[2];
			for i=1:n
				sco_V[i]=0;
			end
			for j=1:G.m
				sco_V[G.u[j]]+=u[G.u[j]]*u[G.v[j]];
				sco_V[G.v[j]]+=u[G.u[j]]*u[G.v[j]];
			end
			for i=1:tot
				x=Int(s[i]);
				sco_V[x]=0;
			end
		else
			d=maxx;
		end
	end
	return s;
end


function greedynode(G,L,k)

	x=argmax(L)[1];
	n=G.n;

	selc2=zeros(k);
	ansp=zeros(k);
	s=union(1:n);
	setdiff!(s,x);
	exist=ones(n);
	orig=zeros(n);
	exist[x]=0;
	xiab=zeros(n);
	t=1;
	for i=1:n
	    if exist[i]>0
	        xiab[i]=t;
	        orig[t]=i;
	        t+=1;
	    end
	end
	selc2[1]=x;
	for i=2:k
	    Lg=L[s,s];
	    e=eigcal(Lg);
	    ansp[i-1]=e[1][1];
	    u=e[2];
	    vect=zeros(n);
	    for t=1:G.m
	        x=G.u[t];
	        y=G.v[t];
	        if x in s && y in s
	            vect[x]+=u[Int(xiab[x])]*u[Int(xiab[y])];
	            vect[y]+=u[Int(xiab[x])]*u[Int(xiab[y])];
	        end
	    end
	    x=argmax(vect);
	    y=x;
	    setdiff!(s,y);
	    exist[Int(y)]=0;
	    selc2[i]=y;
	    xiab=zeros(G.n);
	    t=1;
	    for tt=1:n
	        if exist[tt]>0
	            xiab[tt]=t;
	            orig[t]=tt;
	            t+=1;
	        end
	    end
	end
	return selc2;
end

#=
function getans(L,n,s,k,lset,ls)
	exis = ones(n);
	for i=1:ls
		exis[lset[i]]=0;
	end
	newid=zeros(n);
	newid[1]=exis[1];
	for i=2:n
		newid[i]=newid[i-1]+exis[i];
	end
	orig=zeros(n-ls);
	for i=1:n
		if orig[newid[i]]==0
			orig[newid[i]]=i;
		end
	end
=#
function getans(L,n,s,k)
	D=spzeros(n,n);
	for i=1:k
		x=Int(s[i])
		D[x,x]+=1;
	end

	#F=union(1:n);
	#setdiff!(F,ls);
	#LL=L[F,F];

	e=eigcal(L+D);
	return e[1][1];
end


function exactgreedy(G,L,k,rep)
    x=argmax(L)[1];
    s=zeros(k);
    s[1]=x;
    n=G.n;
    sel=zeros(n);
    sel[x]=1;
    D=spzeros(n,n);
    D[x,x]=1;
    for i=2:k
        maxx=0;
		tmp=0;
        for j=1:G.n
			if sel[j]==0
            	D[j,j]+=1;
				e=eigcal(L+D);
				if e[1][1]>maxx
					maxx=e[1][1]
					tmp=j
				end
				D[j,j]-=1;
			end
        end
        s[i]=tmp;
        sel[tmp]=1;
        D[tmp,tmp]=1;
    end
    return s;
end
