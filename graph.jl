struct Graph
    n :: Int # |V|
    m :: Int # |E|
    u :: Array{Int, 1}
    v :: Array{Int, 1} # uv is an edge
#    w :: Array{Float64, 1} # weight of each edge
    k :: Int
    n0:: Int
    n1:: Int
    s0:: Array{Int, 1}
    s1:: Array{Int, 1}
end

function get_graph(ffname)
    n = 0
    Label = Dict{Int32, Int32}()
    Origin = Dict{Int32, Int32}()
    #E = Set{Tuple{Int32, Int32, Float32}}()

    getID(x :: Int) = haskey(Label, x) ? Label[x] : Label[x] = n += 1

    fname = string("data/",ffname)
    fin = open(fname, "r")
    ##s0,s1
    s0 = Int[]
    s1 = Int[]
    str = readline(fin)
    n0  = parse(Int, str)
    str = readline(fin)
    str = split(str)
    for i = 1 : n0
        push!(s0 , getID(parse(Int, str[i])))
    end
    str = readline(fin)
    n1  = parse(Int, str)
    str = readline(fin)
    str = split(str)
    for i = 1 : n1
        push!(s1 , getID(parse(Int, str[i])))
    end

    str = readline(fin)
    str = split(str)
    #n   = parse(Int, str[1])
    m   = parse(Int, str[3])
    k = parse(Int, str[4])
    u = Int[]
    v = Int[]
    tot = 0
    for i = 1 : m
        str = readline(fin)
        str = split(str)
        x   = parse(Int, str[1])
        y   = parse(Int, str[2])
        if x!=y
            u1 = getID(x)
            v1 = getID(y)
            Origin[u1] = x
            Origin[v1] = y
            push!(u, u1)
            push!(v, v1)
            tot += 1
        end
    end
    for i = 1:k
        push!(u,0)
        push!(v,0)
    end
    close(fin)
    return Graph(n, tot, u, v, k, n0, n1, s0, s1)
end
