
"""An exact MILP model for finding optimal non-overlapping testing regimes."""
function milp_disjoint_model(q, u; T, G=5, verbose=false)
	n = length(q)
	people = 1:n
	M = G * maximum(u) + 1  # Compute big M bound for `implication constraints`

	m = Model(() -> Gurobi.Optimizer(GRB_ENV))
	set_optimizer_attribute(m, "Threads", 8)
	!verbose && set_silent(m)
	# set_optimizer_attribute(m, "TimeLimit", 600)
	# set_optimizer_attribute(m, "Presolve", -1)
    # set_optimizer_attribute(m, "MIPGap", 0.01)
    
	@variable(m, x[1:T, people] >= 0, binary=true)
	@variable(m, s[1:T, 0:n])

	@objective(m, Max, sum(s[t,n] for t in 1:T))

	## Idea: ensure that s[t,i] == u⋅x[t]*prod(q[k]^x[t,k] for k=1:i) for all t and i
	# Encode s[t,0] == x[t]⋅u for all t
	@constraint(m, [t=1:T], s[t,0] == sum(x[t,p] * u[p] for p in people))
    # Encode the implication x[t,k]==0 => s[t,i]==s[t,i-1]
	@constraint(m, [t=1:T, i=1:n], s[t,i] - s[t,i-1] <= M*x[t,people[i]])
    @constraint(m, [t=1:T, i=1:n], s[t,i] - s[t,i-1] >= -M*x[t,people[i]])
	# Encode the implication x[t,i]==1 => s[t,i]==s[t,k-1]*q[k]
	@constraint(m, [t=1:T, i=1:n], s[t,i] - s[t,i-1]*q[i] <= M*(1-x[t,people[i]]))
    @constraint(m, [t=1:T, i=1:n], s[t,i] - s[t,i-1]*q[i] >= -M*(1-x[t,people[i]]))

	# Encode disjoint tests and bounded pool size
	@constraint(m, disjoint[p=people], sum(x[t,p] for t in 1:T) <= 1)  # disjoint tests
	@constraint(m, poolsize[t=1:T], sum(x[t,p] for p in people) <= G)  # bounded pool size

	return m, x
end



"""An exact MILP model for finding optimal overlapping testing regimes."""
function milp_overlap_model(q, u; k, T, G=5, verbose=false)	
	# Compute constants
	n = length(q)  # population size
	people = 1:n
	tests = 1:T
	psets = collect(powerset(tests, 2, k))
	M = G * maximum(u) + 1  # big M bound for the `implication constraints`

	# Create model and set parameters
	m = Model(() -> Gurobi.Optimizer(GRB_ENV))
	set_optimizer_attribute(m, "Threads", 8)
	!verbose && set_silent(m)
	# set_optimizer_attribute(m, "TimeLimit", 600)
	# set_optimizer_attribute(m, "Presolve", -1)
    # set_optimizer_attribute(m, "MIPGap", 0.01)
    
	@variable(m, x[tests, people], Bin)
	@variable(m, con[psets, people], Bin)
	@variable(m, dis[psets, people], Bin)
	@variable(m, sts[tests, 0:n])  # s for tests
	@variable(m, sol[psets, 0:n])  # s for overlaps

	@expression(m, test_contribution, sum(sts[t,n] for t in tests))
	@expression(m, overlap_contribution, sum((-1)^(length(J)+1)*sol[J,n] for J in psets))
	@objective(m, Max, test_contribution + overlap_contribution)

	# Establish relationship between x, con, and dis
	@constraint(m, [J in psets, j in J, i in 1:n], con[J,i] <= x[j,i])
	@constraint(m, [J in psets, i in 1:n], con[J,i] >= 1 - length(J) + sum(x[j,i] for j in J))
	@constraint(m, [J in psets, j in J, i in 1:n], dis[J,i] >= x[j,i])
	@constraint(m, [J in psets, i in 1:n], dis[J,i] <= sum(x[j,i] for j in J))

	## Idea: ensure that sts[t,i] == u⋅x[t]*prod(q[k]^x[t,k] for k=1:i) for all t and i
	# Encode sts[t,0] == x[t]⋅u for all t
	@constraint(m, [t in tests], sts[t,0] == sum(x[t,i] * u[i] for i in people))
    # Encode the implication x[t,k]==0 => sts[t,i]==sts[t,i-1]
	@constraint(m, [t in tests, i in people], sts[t,i] - sts[t,i-1] <= M*x[t,i])
    @constraint(m, [t in tests, i in people], sts[t,i] - sts[t,i-1] >= -M*x[t,i])
	# Encode the implication x[t,i]==1 => sts[t,i]==sts[t,k-1]*q[k]
	@constraint(m, [t in tests, i in people], sts[t,i] - sts[t,i-1]*q[i] <= M*(1-x[t,i]))
    @constraint(m, [t in tests, i in people], sts[t,i] - sts[t,i-1]*q[i] >= -M*(1-x[t,i]))

	## Idea: ensure sol[J,i] == u⋅con[J]*prod(q[k]^dis[J,k] for k=1:i) for all J and i
	@constraint(m, [J in psets], sol[J,0] == sum(con[J,i] * u[i] for i in people))  # ensures sol[J,0] == con[J]⋅u for all t
    # Encode the implication dis[J,k]==0 => sol[J,i]==sol[J,i-1]
	@constraint(m, [J in psets, i in people], sol[J,i] - sol[J,i-1] <= M*dis[J,i])
    @constraint(m, [J in psets, i in people], sol[J,i] - sol[J,i-1] >= -M*dis[J,i])
	# Encode the implication dis[J,i]==1 => sol[J,i]==sol[J,k-1]*q[k]
	@constraint(m, [J in psets, i in people], sol[J,i] - sol[J,i-1]*q[i] <= M*(1-dis[J,i]))
    @constraint(m, [J in psets, i in people], sol[J,i] - sol[J,i-1]*q[i] >= -M*(1-dis[J,i]))

	@constraint(m, disjoint[i in people], sum(x[t,i] for t in tests) <= k)  # k-overlapping tests
	@constraint(m, poolsize[t in tests], sum(x[t,i] for i in people) <= G)  # bounded pool size

	return m, x
end