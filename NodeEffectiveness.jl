include("logw.jl")
include("graph.jl")
include("test.jl")
include("getM.jl")

function exact(G, node_list, pre, w :: IOStream) #G- graph

    n = G.n
    ans = zeros(length(node_list), 1)

    for i in eachindex(node_list)
        ans[i] = test_node(G, node_list[i], pre, w)
    end

    return ans
end

function approx(G, node_list, u, v, w :: IOStream) #G- graph

    n = G.n
    ans = zeros(length(node_list), 1)


    for k in eachindex(node_list)
        edges = Dict{Any, Int}()
        for i in node_list[k]
            for j in G.nbr[i]
                (haskey(edges, (i, j))) && (continue)
                ans[k] += u[i] * v[n+i] + u[n+i] * v[n+j] + u[j] * v[n+j] + u[n+j] * v[n+i]
                edges[(i, j)] = 1
                edges[(j, i)] = 1
            end
        end
    end
    ans /= dot(u, v) * 2

    return ans
end
