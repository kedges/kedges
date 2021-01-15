include("graph.jl")
include("edgecore.jl")
using LinearAlgebra


fname = open("filename.txt", "r")
str   = readline(fname)
n     = parse(Int, str)

for i = 1 : n
    str = readline(fname)
    str = split(str)
    G   = get_graph(str[1])
    #fout = open(string("data/","finnal.ans"),"a")
    #println("Now running file:",str[1],"with eps=",ep);
    for k= 1 : 3
        println("Now running file:",str[1],"with eps=",ep);
        fout = open(string("data/","time.ans"),"a")
        print(fout,str[1]," ");
        ep = 0.1*k;
        greed = gre(G,G.k ,ep,fout);
        println(fout," ");
        close(fout);
    end
    exact=exa(G,G.k,fout);
end

close(fname)
