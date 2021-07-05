include("graph.jl")
include("getM.jl")
include("EigComp.jl")
include("logw.jl")


function test_edge(G :: Graph, edge_list, pre, w :: IOStream) #G- graph

    n = G.n
    ans = 0

    edge = Dict{Tuple{Int, Int}, Int}()
    for (i, j) in edge_list
        edge[(i, j)] = 1
        edge[(j, i)] = 1
    end

    nbr = Array{Array{Int, 1}, 1}()
    for i in 1:n push!(nbr, []) end
    for u in eachindex(G.nbr)
        for v in G.nbr[u]
            if !haskey(edge, (u, v)) push!(nbr[u], v) end
        end
    end

    m = div(sum(length.(nbr)), 2)

    curG = Graph(n, m, nbr)
    curM = getSparseM(curG)
    logw(w, "getM")
    curB = getSparseB(curG)
    #cur = maximum(real.(eigvals(curM)))
    cur, u = eigcomp(curM)
    ans = pre - cur
    return ans
end


function test_node(G :: Graph, node_list, pre, w :: IOStream) #G- graph

    edges = Dict{Any, Int}()
    edge_list = []
    for i in node_list
        for j in G.nbr[i]
            (haskey(edges, (i, j))) && (continue)
            push!(edge_list, (i, j))
            edges[(i, j)] = 1
            edges[(j, i)] = 1
        end
    end
    return test_edge(G, edge_list, pre, w)
end
