include("graph.jl")
include("edgecore.jl")
include("alg1.jl")
#include("alg2.jl")
using LinearAlgebra
using Arpack

fname = open("filename.txt", "r")

str   = readline(fname);
nn     = parse(Int, str);
kmax=5;



for gn=1:nn
    str = readline(fname);
    str = split(str);
    G   = get_graph(str[1]);
    on=G.n;om=G.m;
    Gc=findconnect(G)
    G=Gc;
    n=G.n;
    eps=0.1;
    L=lapsp(G);

    deg=zeros(G.n);
    for i=1:G.n
        deg[i]=L[i,i];
    end
    tt1=time()
    #clo=ClosenessCentrality(G);
    tt2=time()
    #bet=BetweennessCentrality(G);
    tt3=time()
    #A=-L;
    for i=1:n
        #A[i,i]=0;
    end
    #l,r=eigs(A,which=:LM,nev=1);
    #eig=zeros(G.n);
    for i=1:n
    #    eig[i]=abs(r[i]);
    end

    sel_deg=zeros(kmax);
    #sel_clo=zeros(kmax);
    #sel_bet=zeros(kmax);
    #sel_eig=zeros(kmax);



    for i=1:kmax
        x=argmax(deg);
        sel_deg[i]=x;
        deg[x]=0;
#=
        x=argmax(clo);
        sel_clo[i]=x;
        clo[x]=0;

        x=argmax(bet);
        sel_bet[i]=x;
        bet[x]=0;
=#
    #    x=argmax(eig);
#        sel_eig[i]=x;
    #    eig[x]=0;
    end

    #sel_kct=kcent(G,L,kmax);

    t1=time()
    sel_exa=exactgreedy(G,L,kmax,0);
    t2=time()
    sel_gre=greedy(G,L,kmax,0);
    t3=time()

    #sel_ptb=greedy2(G,L,kmax);

    fout = open("buchong.txt", "a")
    println(fout)
    println(fout,str[1],' ',on,' ',om,' ',G.n,' ',G.m," exa");
    println(fout,"exa_time=",t2-t1," fast_time=",t3-t2)
    println(fout,getans(L,n,sel_exa,kmax))
    #println(fout,"clo_time=",tt2-tt1,"bet_time=",tt3-tt2)
    close(fout);


    for k=2:1:0

        n=G.n;
        #s=exactgreedy(G,L,k,0);
        #ans0=getans(L,n,sexa[1:k],k);
        #s=greedy(G,L,k,0);
        #ans1=qpl(k,n,L); #optimum


        #s=greedynode(G,L,k)
        #ans5=getans(L,n,s,k);

        #s=greedy2(G,L,k);
        #ans6=getans(L,n,s,k);
        #for tt=1:1


        #ans1=getans(L,n,sel_deg,k)
        #ans2=0
        ##getans(L,n,sel_exa,k)
        #ans3=getans(L,n,sel_gre,k)
        #ans4=getans(L,n,sel_ptb,k)
        #ans5=getans(L,n,sel_kct,k)
        #ans6=getans(L,n,sel_bet,k)
        #ans7=getans(L,n,sel_clo,k)
        ans8=getans(L,n,sel_eig,k)
        ffout = open("fig4.txt", "a")
        #println(ffout,k,' ',ans1,' ',ans2,' ',ans3,' ',ans4,' ',ans5,' ',ans6,' ',ans7,' ',ans8);
        println(ffout,k,' ',ans8);

        close(ffout)
        #end
    end
end
close(fname);
