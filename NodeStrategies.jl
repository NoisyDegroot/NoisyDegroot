include("logw.jl")
include("graph.jl")
include("test.jl")


function random(G :: Graph, k :: Int, w :: IOStream)
    node_list = []
    for i in 1:k
        u = rand(1:G.n)
        while (u in node_list)
            u = rand(1:G.n)
        end
        push!(node_list, u)
    end
    logw(w, node_list)
    return node_list
end

function topDegree(G :: Graph, k :: Int, w :: IOStream)
    node_list = []
    n = G.n
    for i in 1:n
        Im = length(G.nbr[i])
        push!(node_list, (i, Im))
    end
    sort!(node_list, by = x -> x[2], rev=true)
    node_list = map(t -> t[1], node_list[1:k])
    logw(w, node_list[1:k])
    return node_list
end

function eigenscore(G :: Graph, u, v, k :: Int, w :: IOStream)
    node_list = []
    n = G.n
    for i in 1:n
        Im = u[i] * v[i]
        push!(node_list, (i, Im))
    end
    sort!(node_list, by = x -> x[2], rev=true)
    logw(w, node_list[1:k])
    node_list = map(t -> t[1], node_list[1:k])
    return node_list
end

function closeness(G :: Graph, k :: Int, w :: IOStream)
    node_list = []
    n = G.n
    for i in 1:n
        Im = 0
        vis = fill(-1, G.n)
        q = [i]
        vis[i] = 0
        l, r = 0, 1
        while l < r
        	x = q[l+=1]
        	Im += vis[x]
        	for y in G.nbr[x]
        	  if vis[y] == -1
        	  	vis[y] = vis[x] + 1
        	  	push!(q, y)
        	  	r += 1
        	  end
        	end
        end
        push!(node_list, (i, Im))
    end
    sort!(node_list, by = x -> x[2])
    logw(w, node_list[1:k])
    node_list = map(t -> t[1], node_list[1:k])
    return node_list
end

function betweenness(G :: Graph, k :: Int, w :: IOStream)
    node_list = []
    Im = fill(0.0, G.n)
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
        	  	delta[y] += (sigma[y] / sigma[x]) * (1.0 + delta[x])
        	  end
        	end
        	if x != i
        		Im[x] += delta[x]
        	end
        end
    end
    for i in 1:n
    	push!(node_list, (i, Im[i]))
    end
    sort!(node_list, by = x -> x[2],  rev=true)
    logw(w, node_list[1:k])
    node_list = map(t -> t[1], node_list[1:k])
    return node_list
end

function Imm(G :: Graph, u, v, k :: Int, w :: IOStream)
    node_list = []
    Im = zeros(G.n, 1)
    n = G.n
    for i in 1:n
        for j in G.nbr[i]
            Im[i] += u[i] * v[n+i] + u[n+i] * v[n+j] + u[j] * v[n+j] + u[n+j] * v[n+i]
        end
    end
    for t in 1:k
        i = argmax(Im)[1]
        push!(node_list, i)
        for j in G.nbr[i]
            (Im[j] == -1) && continue
            Im[j] -= u[i] * v[n+i] + u[n+i] * v[n+j] + u[j] * v[n+j] + u[n+j] * v[n+i]
        end
        Im[i] = -1
    end
    logw(w, node_list)
    return node_list
end


function optimal(G :: Graph, k :: Int, pre_max_u, w :: IOStream)
    node_list, ans = [], 0
    cur_node_list = []

    function dfs(k :: Int, cur :: Int)
        if k == 0
            cur_ans = test_node(G, cur_node_list, pre_max_u, w)
            #logw(w, "node: ", cur_node_list, "\tans: ", cur_ans)
            if cur_ans > ans
                ans = cur_ans
                node_list = copy(cur_node_list)
            end
            return
        end
        for i in cur:G.n
            push!(cur_node_list, i)
            dfs(k-1, i+1)
            pop!(cur_node_list)
        end
    end

    dfs(k, 1)
    logw(w, node_list)
    return node_list, ans
end
