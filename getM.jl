using LinearAlgebra
using SparseArrays
include("graph.jl")

function getA(G :: Graph)
    A :: Matrix{Float64} = zeros(G.n, G.n)
    for i in eachindex(G.nbr)
        A[CartesianIndex.(tuple.(i, G.nbr[i]))] .=  1.0
    end
    return A
end

function getM(G :: Graph)
    A :: Matrix{Float64} = zeros(G.n, G.n)
    D :: Matrix{Float64} = Diagonal(size.(G.nbr, 1))
    O :: Matrix{Float64} = zeros(G.n, G.n)
    II :: Matrix{Float64} = Matrix{Float64}(I, G.n, G.n)
    for i in eachindex(G.nbr)
        A[CartesianIndex.(tuple.(i, G.nbr[i]))] .=  1.0
    end
    M :: Matrix{Float64} = vcat(hcat(O, D - II), hcat(-II, A))
    return M
end

function getSparseA(G :: Graph)
    I, J, V :: Array{Float64, 1} = [], [], []

    for i in eachindex(G.nbr)
        for j in G.nbr[i]
            push!(I, i)
            push!(J, j)
	    push!(V, 1)
        end
    end

    return sparse(I, J, V, G.n, G.n)
end

function getSparseM(G :: Graph)
    I, J, V :: Array{Float64, 1} = [], [], []

    for i in eachindex(G.nbr)
        for j in G.nbr[i]
            push!(I, G.n + i)
            push!(J, G.n + j)
            push!(V, 1)
        end
        push!(I, i)
        push!(J, G.n + i)
        push!(V, length(G.nbr[i]) - 1)
        push!(I, G.n + i)
        push!(J, i)
        push!(V, -1)
    end

    return sparse(I, J, V, 2 * G.n, 2 * G.n)
end

function getSparseB(G :: Graph)
    I, J, V :: Array{Float64, 1} = [], [], []
    edge_list = Dict{Tuple{Int, Int}, Int}()
    edgelist = []
    for i in eachindex(G.nbr)
        for j in G.nbr[i]
            push!(edgelist, (i, j))
        end
    end
    for i in eachindex(edgelist)
        edge_list[edgelist[i]] = i
    end
    for (i, j) in edgelist
        for k in G.nbr[j]
            (k == i) && continue
            push!(I, edge_list[(i, j)])
            push!(J, edge_list[(j, k)])
            push!(V, 1)
        end
    end

    return sparse(I, J, V, 2 * G.m, 2 * G.m), edgelist
end
