include("graph.jl")
include("logw.jl")
include("getM.jl")
include("EigComp.jl")
include("test.jl")


#using StatsBase


datadir = string(ARGS[1], "/")
outFName=string(ARGS[1], ".txt")
w = open(outFName, "w")

if ARGS[2] == "node"
	include("NodeStrategies.jl")
	include("NodeEffectiveness.jl")
elseif ARGS[2] == "edge"
	include("Effectiveness.jl")
	include("Strategies.jl")
end


for rFile in filter(x -> endswith(x, ".txt"), readdir(string(datadir)))
    logw(w, "reading graph from edges list file ", rFile)
    G = read_file(string(datadir, rFile))
    logw(w, "finished reading graph")
    logw(w, "LCC.n: ", G.n, "\t LCC.m: ", G.m)

	logw(w, "")
    logw(w, "start")
    M = getSparseM(G)
	logw(w, "get M")
    logw(w, "computing the max eigenvalue of M")

	start_time = time()
	#logw(w, "eigvals: ", eigvals(M))
	#logw(w, "leading eigval for B: ", getSparseB(G)[1])
	#pre_max_v = maximum(real.(eigvals(M)))
	pre_max_u, u = eigcomp(M')
	logw(w, "max left eig:", pre_max_u)
	pre_max_v, v = eigcomp(M)
	logw(w, "max right eig:", pre_max_v)
	A = getSparseA(G)
	pre_max_uA, uA = eigcomp(A)
	logw(w, "max right eig for A:", pre_max_uA)
	elapsed_time = time() - start_time

	logw(w, "")
#=
	# Effectiveness
	d = length.(G.nbr)
	if ARGS[2] == "edge"
		edges = []
		W = zeros(G.m, 1)
		cur = 0
		for i in 1:G.n
			for j in G.nbr[i]
				(i > j) && continue
				push!(edges, (i, j))
				cur += 1
				W[cur] = d[i] * d[j]
			end
		end
	elseif ARGS[2] == "node"
		nodes = []
		mx = maximum(d)
		W = zeros(G.n, 1)
		for i in 1:G.n
			W[i] = 1 / (mx - d[i] + 1)
			push!(nodes, i)
			#(d[i] > 3) ? (W[i] = 0) : (W[i] = 5-W[i]+1)
		end
	end

	cnt = 1
	for k in 2:2
		list = []
		for i in 1:cnt
			if ARGS[2] == "edge"
				edge = StatsBase.sample(edges, ProbabilityWeights(W[:]), 1; replace=false, ordered=true)
				push!(list, edge)
			elseif ARGS[2] == "node"
				node = StatsBase.sample(nodes, ProbabilityWeights(W[:]), k; replace=false, ordered=true)
				println((node, d[node]))
				push!(list, node)
			end
		end
		logw(w, "")
		logw(w, "list: ", list)
		logw(w, "length of list: ", cnt, "\t k: ", k)

		exact_rst = exact(G, list, pre_max_u, w)
		logw(w, "cent: ", exact_rst)

	  	approx_rst = approx(G, list, u, v, w)
		logw(w, "cent: ", approx_rst)

		mRelStdEr, maxRelStdEr, tot = 0, 0, 0
		for i in 1:cnt
			if abs(exact_rst[i]) < 1e-7 continue end
			tot += 1
			mRelStdEr += abs(exact_rst[i] - approx_rst[i])/exact_rst[i]
			maxRelStdEr = max(maxRelStdEr, abs(exact_rst[i] - approx_rst[i])/exact_rst[i])
		end
		logw(w, "valid_cnt: ", tot)
	    logw(w, "approx_mRelStdEr:", mRelStdEr/tot)
	    logw(w, "approx_maxRelStdEr:", maxRelStdEr)
	end
=#

	# Performances

	if ARGS[2] == "node"
		test = test_node
	elseif ARGS[2] == "edge"
		test = test_edge
	end

	k_max = 20
#=
	# ExactImm
	ExactImm_list = ExactImm(G, pre_max_u, k_max, w)
=#
	# Imm
	Imm_list = Imm(G, u, v, k_max, w)

	# ImmA
	ImmA_list = ImmA(G, uA, k_max, w)
#=	
	# Random
	random_list = random(G, k_max, w)	

	# Eigen-score
	eigsc_list = eigenscore(G, u[G.n+1:end], v[G.n+1:end], k_max, w)
	
	# Betweenness
	betweenness_list = betweenness(G, k_max, w)
			
	if ARGS[2] == "edge"
		# Top-degree-product
		topDegreeProduct_list = topDegreeProduct(G, k_max, w)

		# Top-degree-sum
		topDegreeSum_list = topDegreeSum(G, k_max, w)
	end

	if ARGS[2] == "node"
		# Top-degree
		topDegree_list = topDegree(G, k_max, w)
		
		# Closeness
		closeness_list = closeness(G, k_max, w)
	end
=#
	for k in 1:20
		logw(w, "k: ", k)
#=	
		ExactImm_ans = test(G, ExactImm_list[1:k], pre_max_u, w)
		logw(w, "Eigen-drop for ExactImm: ", ExactImm_ans)
		logw(w, "")
=#
		Imm_ans = test(G, Imm_list[1:k], pre_max_u, w)
		logw(w, "Eigen-drop for Imm: ", Imm_ans)
		logw(w, "")

		ImmA_ans = test(G, ImmA_list[1:k], pre_max_u, w)
		logw(w, "Eigen-drop for ImmA: ", ImmA_ans)
		logw(w, "")
#=		
		random_ans = test(G, random_list[1:k], pre_max_u, w)
		logw(w, "Eigen-drop for Random: ", random_ans)
		logw(w, "")
		
		eigsc_ans = test(G, eigsc_list[1:k], pre_max_u, w)
		logw(w, "Eigen-drop for Eigen-score: ", eigsc_ans)
		logw(w, "")
	
		betweenness_ans = test(G, betweenness_list[1:k], pre_max_u, w)
		logw(w, "Eigen-drop for Betweenness: ", betweenness_ans)
		logw(w, "")
			
		if ARGS[2] == "edge"
			topDegreeProduct_ans = test(G, topDegreeProduct_list[1:k], pre_max_u, w)
			logw(w, "Eigen-drop for TopDegreeProduct: ", topDegreeProduct_ans)
			logw(w, "")

			topDegreeSum_ans = test(G, topDegreeSum_list[1:k], pre_max_u, w)
			logw(w, "Eigen-drop for TopDegreeSum: ", topDegreeSum_ans)
			logw(w, "")
		end
=#
		if ARGS[2] == "node"
			topDegree_ans = test(G, topDegree_list[1:k], pre_max_u, w)
			logw(w, "Eigen-drop for TopDegree: ", topDegree_ans)
			logw(w, "")
			
			closeness_ans = test(G, closeness_list[1:k], pre_max_u, w)
			logw(w, "Eigen-drop for Closeness: ", closeness_ans)
			logw(w, "")		
		end
#=
		# Optimal
		list, optimal_ans = optimal(G, k, pre_max_u, w)
		logw(w, "Eigen-drop for Optimal: ", optimal_ans)
		logw(w, "")
=#
		logw(w, "")
    	logw(w, String(fill('*', 60)))
    	logw(w, "")
	end

    logw(w, "")
    logw(w, String(fill('*', 80)))
    logw(w, "")
end
close(w)
