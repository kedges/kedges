include("graph.jl")
include("edgecore.jl")
using LinearAlgebra
#using BenchmarkTools


fname = open("filename.txt", "r")
str   = readline(fname)
n     = parse(Int, str)

eps=0.1

for i = 1:n
    str = readline(fname)
    str = split(str)
    G   = get_graph(str[1])
    #fout = open(string("data/","finnal.ans"),"a")
    #精确解和近似解 k=1 - G.k
    for hh=1:3
    eps = 0.1*hh;
    println("Now running file:",str[1],"with eps=",eps);
    for k= G.k/2 : G.k
        fout = open(string("data/","finnal.ans"),"a")
        T0=time();
        exact = exa(G,k);
        T1=time();
        #@benchmark
        greed = gre(G,k ,eps);
        T2=time();
        println(fout,str[1]," ,k=",k," exact=",exact," ,greedy =", greed, " ,eps=",eps," ,error=",1-greed/exact,", exact time=",T1-T0,", greedy time=",T2-T1);
        close(fout);
    end
    #close(fout)
    end

end

close(fname)
