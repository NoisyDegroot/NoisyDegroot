include("logw.jl")
include("graph.jl")
include("test.jl")


function random(G :: Graph, k :: Int, w :: IOStream)
    edge_list = []
    for i in 1:k
        (u, v) = rand_edge(G)
        while ((u, v) in edge_list) || ((v, u) in edge_list)
            (u, v) = rand_edge(G)
        end
        push!(edge_list, (u, v))
    end
    logw(w, edge_list)
    return edge_list
end

function topDegreeProduct(G :: Graph, k :: Int, w :: IOStream)
    edge_list = []
    n = G.n
    for i in 1:n
        for j in G.nbr[i]
            (i > j) && (continue)
            Im = length(G.nbr[i]) * length(G.nbr[j])
            push!(edge_list, ((i, j), Im))
        end
    end
    sort!(edge_list, by = x -> x[2], rev=true)
    logw(w, edge_list[1:k])
    edge_list = map(t -> t[1], edge_list[1:k])
    return edge_list
end

function topDegreeSum(G :: Graph, k :: Int, w :: IOStream)
    edge_list = []
    n = G.n
    for i in 1:n
        for j in G.nbr[i]
            (i > j) && (continue)
            Im = length(G.nbr[i]) + length(G.nbr[j])
            push!(edge_list, ((i, j), Im))
        end
    end
    sort!(edge_list, by = x -> x[2], rev=true)
    logw(w, edge_list[1:k])
    edge_list = map(t -> t[1], edge_list[1:k])
    return edge_list
end

function eigenscore(G :: Graph, u, v, k :: Int, w :: IOStream)
    edge_list = []
    n = G.n
    for i in 1:n
        for j in G.nbr[i]
            (i > j) && (continue)
            Im = u[i] * v[j] + u[j] * v[i]
            push!(edge_list, ((i, j), Im))
        end
    end
    sort!(edge_list, by = x -> x[2], rev=true)
    logw(w, edge_list[1:k])
    edge_list = map(t -> t[1], edge_list[1:k])
    return edge_list
end

function betweenness(G :: Graph, k :: Int, w :: IOStream)
    edge_list = []
    Im = Dict{Tuple{Int, Int}, Float64}()
    n = G.n
    for i in 1:n
        vis = fill(-1, G.n)
        sigma = fill(1.0, G.n)
        delta = fill(0.0, G.n)
        q = [i]
        vis[i] = 0
        l, r = 0, 1
        while l < r
        	x = q[l+=1]
        	for y in G.nbr[x]
        	  if vis[y] == -1
        	  	vis[y] = vis[x] + 1
        	  	push!(q, y)
        	  	r += 1
        	  elseif vis[y] == vis[x] + 1
        	  	sigma[y] += sigma[x]
        	  end
        	end
        end
        for j in n:-1:1
        	x = q[j]
        	for y in G.nbr[x]
        	  if vis[x] == vis[y] + 1
        	  	c = (sigma[y] / sigma[x]) * (1.0 + delta[x])
        	  	delta[y] += c
        	  	(haskey(Im, (min(x, y), max(x, y)))) ? (Im[(min(x, y), max(x, y))] += c) : (Im[(min(x, y), max(x, y))] = c)
        	  end
        	end
        end
    end
    for (i, j) in keys(Im)
    	push!(edge_list, ((i, j), Im[(i, j)]))
    end
    sort!(edge_list, by = x -> x[2],  rev=true)
    logw(w, edge_list[1:k])
    edge_list = map(t -> t[1], edge_list[1:k])
    return edge_list
end

function Imm(G :: Graph, u, v, k :: Int, w :: IOStream)
    edge_list = []
    n = G.n
    for i in 1:n
        for j in G.nbr[i]
            (i > j) && (continue)
            Im = u[i] * v[n+i] + u[n+i] * v[n+j] + u[j] * v[n+j] + u[n+j] * v[n+i]
            push!(edge_list, ((i, j), Im))
        end
    end
    sort!(edge_list, by = x -> x[2], rev=true)
    # logw(w, edge_list[1:k])
    edge_list = map(t -> t[1], edge_list[1:k])
    return edge_list
end

function ImmA(G :: Graph, uA, k :: Int, w :: IOStream)
    edge_list = []
    n = G.n
    for i in 1:n
        for j in G.nbr[i]
            (i > j) && (continue)
            Im = uA[i] * uA[j]
            push!(edge_list, ((i, j), Im))
        end
    end
    sort!(edge_list, by = x -> x[2], rev=true)
    # logw(w, edge_list[1:k])
    edge_list = map(t -> t[1], edge_list[1:k])
    return edge_list
end

function ExactImm(G :: Graph, pre, k :: Int, w :: IOStream)
    edge_list = []
    n = G.n
    for i in 1:n
        for j in G.nbr[i]
            (i > j) && (continue)
            Im = test_edge(G, [(i, j)], pre, w)
            push!(edge_list, ((i, j), Im))
        end
    end
    sort!(edge_list, by = x -> x[2], rev=true)
    logw(w, edge_list[1:k])
    edge_list = map(t -> t[1], edge_list[1:k])
    return edge_list
end

function optimal(G :: Graph, k :: Int, pre_max_u, w :: IOStream)
    edge_list, ans = [], 0
    cur_edge_list = []
    edges = []
    for i in 1:G.n
        for j in G.nbr[i]
            (i < j) && (push!(edges, (i, j)))
        end
    end

    function dfs(k :: Int, cur :: Int)
        if k == 0
            cur_ans, elapsed_time = test_edge(G, cur_edge_list, pre_max_u, w)
            logw(w, "edge: ", cur_edge_list, "\tans: ", cur_ans)
            if cur_ans > ans
                ans = cur_ans
                edge_list = copy(cur_edge_list)
            end
            return
        end
        for i in cur:G.m
            push!(cur_edge_list, edges[i])
            dfs(k-1, i+1)
            pop!(cur_edge_list)
        end
    end

    dfs(k, 1)
    logw(w, edge_list)
    return edge_list, ans
end
