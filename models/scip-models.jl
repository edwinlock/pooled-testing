using JuMP, SCIP, AmplNLWriter, Bonmin_jll
using Combinatorics

"""
# Instructions
1. Install Julia (from www.julialang.org)
2. Install Julia packages: start Julia in the terminal, press ] to enter package mode and type `add JuMP, SCIP`.
3. Run this script by typing `julia testing-scip.jl` in the terminal.
"""


"""
Build cluster-based model for SCIP to solve.
Inputs are three vectors indexed by cluster, as well as number of tests T
and pool size bound G.

u::Vector - utility for each cluster
q::Vector - avg. probability of being healthy for each cluster
n::Vector - size of each cluster
T::Int 	  - number of tests
G::Int    - pooled test size
"""
function cluster_scip(u, q, n; T=1, G=10)
	@assert length(u) == length(q) == length(n) # check input is consistent
	C = length(n)  # number of clusters C

	# Create model and set parameters
	# m = Model(SCIP.Optimizer)
	m = Model(() -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe))
	# JuMP.set_optimizer_attribute(m, "display/verblevel", 0)
	# set_optimizer_attribute(m, "presolving/maxrounds", 0)

	# Define variables
	@variable(m, x[1:T, 1:C] >= 0, Int)  # constraint (2h)
	@variable(m, l[1:T])
	@variable(m, w[1:T])
	@variable(m, y[1:T] >= 0)
	@variable(m, z[1:T] >= 0)

	# Set objective
	@objective(m, Max, sum(w[t] for t in 1:T))

	# Add constraints
	@NLconstraint(m, [t=1:T], w[t] <= exp(l[t]))  # (2b)
	@constraint(m, [t=1:T], l[t] == y[t] + sum(x[t,i] * log(q[i]) for i in 1:C))  # (2c)
	@NLconstraint(m, [t=1:T], y[t] <= log(z[t]))  # (2d)
	@constraint(m, [t=1:T], z[t] == sum(x[t,i]*u[i] for i in 1:C))  # (2e)
	@constraint(m, [i=1:C], sum(x[t,i] for t in 1:T) <= n[i])  # (2f)
	@constraint(m, [t=1:T], 1<= sum(x[t,i] for i in 1:C) <= G)  # (2g)

	# Return model and variables
	return m, x
end

disjoint_scip(u, q; T=1, G=10) = cluster_scip(u, q, ones(length(u)); T=T, G=G)


"""
We now allow people to lie in two tests but no longer have clusters.
"""
function overlaps_scip(u, q; T, G)
	@assert length(u) == length(q)  # check input is consistent
	n = length(u) # size of population
	testpairs = collect(combinations(1:T, 2))  # compute all possible pairs of tests
    TP = length(testpairs)

	# Create model and set parameters
	m = Model(SCIP.Optimizer)
	# m = Model(() -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe))
	# JuMP.set_optimizer_attribute(m, "display/verblevel", 0)
	# set_optimizer_attribute(m, "presolving/maxrounds", 0)

	# Single test variables
	@variable(m, x[1:T, 1:n], binary=true)
	@variable(m, l[1:T])
	@variable(m, w[1:T])
	@variable(m, y[1:T])
	@variable(m, z[1:T])

	# Test pair variables
	@variable(m, a[1:TP, 1:n]>=0, binary=true)
    @variable(m, b[1:TP, 1:n]>=0, binary=true)
	@variable(m, wp[1:TP])
    @variable(m, lp[1:TP])
    @variable(m, yp[1:TP])
    @variable(m, zp[1:TP])
    
	# Set objective
	@objective(m, Max, sum(w[t] for t in 1:T) + sum(-wp[p] for p in 1:TP))

    # Single test constraints
    @NLconstraint(m, [t=1:T], w[t] == exp(l[t]))
	@constraint(m, [t=1:T], l[t] == y[t] + sum(x[t,i] * log(q[i]) for i in 1:n))
	@NLconstraint(m, [t=1:T], y[t] == log(z[t]))
	@constraint(m, [t=1:T], z[t] == sum(x[t,i]*u[i] for i in 1:n))

	# Test pair constraints
    @NLconstraint(m, [p=1:TP], wp[p] == exp(lp[p]))
	@constraint(m, [p=1:TP], lp[p] == yp[p] + sum(a[p,i] * log(q[i]) for i in 1:n))
	@NLconstraint(m, [p=1:TP], yp[p] == log(zp[p]))
	@constraint(m, [p=1:TP], zp[p] == sum(b[p,i]*u[i] for i in 1:n))

	# Constraints governing values of a and b
	for p in 1:TP
		s, t = testpairs[p]
		println(p, " ", s, " ", t)
		# Ensure that a[p,i] == (x[s,i] OR x[t,i])
		@constraint(m, [i=1:n], x[s,i]+x[t,i]<=2*a[p,i])
		@constraint(m, [i=1:n], x[s,i]+x[t,i]>=a[p,i])
		# Ensure that b[p,i] == (x[s,i] AND x[t,i])
		@constraint(m, [i=1:n], x[s,i]+x[t,i]<=1+b[p,i])
		@constraint(m, [i=1:n], x[s,i]+x[t,i]>=2b[p,i])
	end

    # Pool size and 2-overlap constraints
	@constraint(m, [i=1:n], sum(x[t,i] for t in 1:T) <= 2)  # 2-overlaps
	@constraint(m, [t=1:T], sum(x[t,i] for i in 1:n) <= G)  # pool size

	# Return model and variables
	return m, x, a, b
end



function scip_overlap_model(q, u; k, T, G)
	# TODO: add asserts to check validity of inputs

	n = length(q)  # population size
	psets = collect(powerset(1:T, 1, k))
	ε = 1e-10
	m = Model(SCIP.Optimizer)

	@variable(m, x[1:T, 1:n], Bin)
	@variable(m, con[psets, 1:n], Bin)
	@variable(m, dis[psets, 1:n], Bin)
	@variable(m, w[psets] >= 0)
	@variable(m, l[psets])
	@variable(m, y[psets])
	@variable(m, z[psets])
	@variable(m, b[psets], Bin)

	@objective(m, Max, sum((-1)^(length(J)+1)*w[J] for J in psets))

	# Establish relationship between x, con, and dis
	@constraint(m, [J in psets, j in J, i in 1:n], con[J,i] <= x[j,i])
	@constraint(m, [J in psets, i in 1:n], con[J,i] >= 1 - length(J) + sum(x[j,i] for j in J))
	@constraint(m, [J in psets, j in J, i in 1:n], dis[J,i] >= x[j,i])
	@constraint(m, [J in psets, i in 1:n], dis[J,i] <= sum(x[j,i] for j in J))

	@NLconstraint(m, [J in psets], w[J] == exp(l[J]))
	@constraint(m, [J in psets], l[J] == y[J] + sum(dis[J,i] * log(q[i]) for i in 1:n))
	@constraint(m, [J in psets], z[J] == sum(con[J,i]*u[i] for i in 1:n))
	@NLconstraint(m, [J in psets], y[J] == log(z[J]))
	# @NLconstraint(m, [J in psets], y[J] == -1e6*b[J] + log(z[J])*(1-b[J]))
	# @constraint(m, [J in psets], b[J]*z[J] <= ε)
	# @constraint(m, [J in psets], z[J] >= (1-b[J]) * ε)

	@constraint(m, [i in 1:n], sum(x[t,i] for t in 1:T) <= k)  # k-overlaps
	@constraint(m, [t in 1:T], sum(x[t,i] for i in 1:n) <= G)  # pool size

	return m, x, dis, con, w, l, y, z, b
end



# # pop = Dict(1 => (0.5, 1) , 2 => (0.5, 1), 3 => (1., 1)
# n = 2
# # pop = generate_instance(n, 0.1:0.1:1, 1:10)
# pop = Dict(1 => (0.5, 1) , 2 => (0.5, 1))
# q, u = pop2vec(pop)
# T = 2
# G = 1

# m, x = scip_overlap_model(q, u; k=2, T=T, G=G)
# optimize!(m)
# value.(x)
# w, pools = retrieve(m, x, T, n)
# welfare(pools, pop)
# m, x, dis, con, w, l, y, z, b = scip_overlap_model(q, u; k=2, T=T, G=G)
# point = Dict(
# 	x[1,1] => 1.,
# 	x[1,2] => 0.,
# 	x[2,1] => 0.,
# 	x[2,2] => 1.,
# 	b[[1]] => 0.,
# 	b[[2]] => 0.,
# 	b[[1,2]] => 1.,
# 	con[[1],1] => 1.,
# 	con[[1],2] => 0.,
# 	con[[2],1] => 0.,
# 	con[[2],2] => 1.,
# 	con[[1,2],1] => 0.,
# 	con[[1,2],2] => 0.,
# 	dis[[1],1] => 1.,
# 	dis[[1],2] => 0.,
# 	dis[[2],1] => 0.,
# 	dis[[2],2] => 1.,
# 	dis[[1,2],1] => 1.,
# 	dis[[1,2],2] => 1.,
# 	w[[1]] => 0.5,
# 	w[[2]] => 0.5,
# 	w[[1,2]] => 0.,
# 	l[[1]] => log(0.5),
# 	l[[2]] => log(0.5),
# 	l[[1,2]] => -Inf,
# 	y[[1]] => log(1),
# 	y[[2]] => log(1),
# 	y[[1,2]] => -Inf,
# 	z[[1]] => 1.,
# 	z[[2]] => 1.,
# 	z[[1,2]] => 0.
# )
# primal_feasibility_report(m, point)