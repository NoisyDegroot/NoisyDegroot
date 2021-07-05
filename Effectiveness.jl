include("logw.jl")
include("graph.jl")
include("test.jl")
include("getM.jl")

function exact(G, edge_list, pre, w :: IOStream) #G- graph

    n = G.n
    ans = zeros(length(edge_list), 1)
    for i in eachindex(edge_list)
        ans[i] = test_edge(G, edge_list[i], pre, w)
    end
    return ans
end

function approx(G, edge_list, u, v, w :: IOStream) #G- graph

    n = G.n
    ans = zeros(length(edge_list), 1)
    for k in eachindex(edge_list)
        for (i, j) in edge_list[k]
            ans[k] += u[i] * v[n+i] + u[n+i] * v[n+j] + u[j] * v[n+j] + u[n+j] * v[n+i]
        end
    end
    ans /= dot(u, v)
    return ans
end
