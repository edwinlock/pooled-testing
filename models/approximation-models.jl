using Roots, Combinatorics

"""
Build cluster-based MILP model for Gurobi to solve. Finds an approximately
optimal testing allocation.
Inputs are three vectors indexed by cluster, as well as number of tests T
and pool size bound G, and accuracy parameter K.

# Arguments
- `u::Vector`: utility for each cluster
- `q::Vector`: avg. probability of being healthy for each cluster
- `n::Vector`: size of each cluster
- `T::Int`: number of tests
- `G::Int=5`: pooled test size
- `K::Int=15`: number of segments of piecewise-linear fn approximating exp constraint

"""
function approx_disjoint_model(q, u, n; T, G=5, K=15, verbose=false)
	# Verify that input is consistent
	@assert length(u) == length(q) == length(n) "Input vectors have different lengths."
	@assert K >= 1  "Number of segments for approximating exp must be at least 1."
	@assert T <= sum(n)  "Number of tests cannot exceed number of people in population."
	@assert all(isinteger.(u)) && all(u .>= 0)  "Utilities must be strictly positive integers."
	@assert all(0 .<= q .<= 1)  "Probabilities must lie strictly between 0 and 1."

	# Compute some constants
	C = length(n)  # number of clusters C
	L, U = minimum(u), G*maximum(u)  # lower and upper bound for z[t] = x[t]⋅u
	# Lower and upper bound for l[t] = log(x[t]⋅u) + sum(x[t,i]*log(q[i])
	A, B = minimum(log.(u)) + G*minimum(log.(q)), log(G*maximum(u)) + maximum(log.(q))

	# Create model and set parameters
	m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV)))
	set_optimizer_attribute(m, "Threads", 8)
	# "Presolve" => 0, "OutputFlag" => 0, "MIPGap" => 0.01, "TimeLimit" => 600
	!verbose && set_silent(m)

	# Define variables
	@variable(m, x[1:T, 1:C] >= 0, Int)
	@variable(m, w[1:T])
	@variable(m, l[1:T])
	@variable(m, y[1:T])
	@variable(m, z[1:T])

	# Set objective
	@objective(m, Max, sum(w[t] for t in 1:T))

	# Add constraints
	@constraint(m, [i=1:C], sum(x[t,i] for t in 1:T) <= n[i])  # tests must be disjoint
	@constraint(m, [t=1:T], 1 <= sum(x[t,i] for i in 1:C) <= G)  # pool size between 1 and G

	# Deal with log welfare constraints: l[t] == log(u ̇x[t]) + x[t] ̇log.(q)
	@constraint(m, [t=1:T], l[t] == y[t] + sum(x[t,i] * log(q[i]) for i in 1:C))  # log welfare
	# Constraints to ensure y[t] <= log(z[t])
	@constraint(m, [t=1:T], z[t] == sum(x[t,i] * u[i] for i in 1:C))  # utility sums z[t] = x[t]⋅u
	# Use indictator variables to capture value of z[t]:
	# z[t] is an integer in [L, U], so let zind[t,k] = 1 if z[t] = k and 0 otherwise.
	@variable(m, zind[1:T, L:U], binary=true)
	@constraint(m, [t=1:T], 1 == sum(zind[t,k] for k in L:U))  # exactly one zind entry is 1
	@constraint(m, [t=1:T], z[t] == sum(k*zind[t,k] for k in L:U))
	@constraint(m, [t=1:T], y[t] == sum(log(k)*zind[t,k] for k in L:U))

	# Deal with exp constraints w[t] = exp(l[t])
	# Variables for approximating exp constraint with a piecewise-linear function
	@variable(m, lind[1:T, 1:K], binary=true)
	@variable(m, v[1:T, 1:K])
	# Approximate w[t] <= exp(l[t]) using piecewise-linear function f with K segments on domain [A,B]
	a, b, c, ε = upper_bounds(A, B, K)  # compute optimal segmentation of interval [A, B]
	# Use indicator variables `lind[t,k]` to capture in which segment the value of l[t] lies
	# and let v[t,k] = l[t] if l[t] lies in (c[k], c[k+1]) and v[t,k] = 0 otherwise.
	@constraint(m, [t in 1:T], 1 == sum(lind[t,k] for k in 1:K))  # exactly one lind entry is 1
	@constraint(m, [t in 1:T, k=1:K], c[k]*lind[t,k] <= v[t,k])
	@constraint(m, [t in 1:T, k=1:K], v[t,k] <= c[k+1]*lind[t,k])
	@constraint(m, [t in 1:T], l[t] == sum(v[t,k] for k in 1:K))  # l[t] = Σv[t,k]
	# Ensure that w[t] <= f(l[t])
	@constraint(m, [t in 1:T], w[t] <= sum(a[k]*v[t,k] + b[k]*lind[t,k] for k in 1:K))

	# Return model
	return m, x, T*ε
end


"""
Build individual-based MILP model for Gurobi to solve. Finds an approximately
optimal testing allocation. Inputs are three vectors indexed by individual, as
well as number of tests T and pool size bound G, and accuracy parameter K.

# Arguments
- `u::Vector`: utility for each individual
- `q::Vector`: avg. probability of being healthy for each individual
- `T::Int`: number of tests
- `G::Int=5`: pooled test size
- `K::Int=15`: number of segments of piecewise-linear fn approximating exp constraint

"""
approx_disjoint_model(q, u; T, G=5, K=15, verbose=false) = approx_disjoint_model(q, u, ones(length(u)); T=T, G=G, K=K, verbose=verbose)


"""
Compute maximum difference between segment (l, exp(l)) to (r, exp(r))
and exp(x) on the interval [l,r].
"""
function Δ(l, r)
	r <= l && return 0.0
	a = (exp(l) - exp(r)) / (l - r)
	a == 0 && return 0.0  # happens if l and r are sufficiently similar
	b = exp(r) - a*r
	result = a*log(a) + b - a   # maximum difference, derived from first order conditions
	return max(0, result)  # slight hack to avoid numerical inaccuracies
end


"""
Build a partition of K segments starting from A such that the
first segment is [A, A+r1] and all segments have the same error
ε identical to the error of the first segment.
"""
function partition(A, K, r1)
	@assert r1 >= 0
	@assert K > 0
	c = fill(float(A), K+1)
	r1 ≈ 0 && return c  # nothing to do if first segment has size 0
	ε = Δ(A, A + r1)  # error of the first segment [A, A+r1]
	for k in 1:K
		l = c[k]
		# To define the bracket for the root finder, we make the reasonable
        # assumption that the interval will be no larger than r1. (This can
        # be proved easily, I believe.)
		r = fzero(x->Δ(l,x)-ε, l, l+r1+1)  # Find r such that Δ(l,r) = ε.
		c[k+1] = r
	end
	return c
end


"""
Find the optimal partition of [A, B] into K segments. Proceeds by searching
for the right size for the first segment: the size `r1` is right when
`partition(A, K, r1)` ends (approximately) at `B`.
"""
function upper_bounds(A, B, K)
	@assert A <= B
	if A ≈ B
		a = fill(exp(A), K)
		c = fill(float(A), K+1)
		b = a .* (1 .- c[1:K])
		return a, b, c
	end
	r1 = fzero(r->partition(A, K, r)[end]-B, 0, B-A+1)
	c = partition(A, K, r1)
	c[K+1] = B  # to clean things up a bit
	# Compute linearisation
	a, b = zeros(K), zeros(K)
	for k in 1:K
		a[k] = (exp(c[k+1]) - exp(c[k])) / (c[k+1] - c[k])  # slope
		b[k] = exp(c[k+1]) - a[k]*c[k+1]  # residual
	end
	error = maximum(a.*log.(a).+b.-a)
	return a, b, c, error
end


function lower_bounds(A, B, K)
	@assert A <= B
	if A ≈ B
		a = fill(exp(A), K+1)
		c = fill(float(A), K+1)
		b =  a .* (1 .- c)
		return a, b, c
	end
	error = fzero(x->tangents(A, K, x)[3][end]-B, 1)
	a, b, c = tangents(A, K, error)
	c[K+1] = B  # to clean things up a bit
	return a, b, c, error
end


function tangents(A, K, ε)
	@assert K >= 1 "number of segments must be at least 1"
	c = fill(float(A), K+1)
	a = fill(exp(A), K+1)
	b =  a .* (1 .- c)
	ε < 1e-10 && return a, b, c  # nothing to do if first segment has size 0
	for i in 2:K+1
		x = fzero(x->exp(x)-a[i-1]*x-b[i-1]-ε, c[i-1]+1)
		c[i] = fzero(c->a[i-1]*x+b[i-1]-exp(c)*(x+1-c), x+1)
		a[i] = exp(c[i])
		b[i] = a[i]*(1-c[i])
	end
	return a, b, c
end


"""
Unfinished implementation for finding k-overlapping testing regimes.
"""
function approx_model(q, u; k=1, T, G=5, K=15, verbose=false)
	# Verify that input is consistent
	@assert length(u) == length(q) "Input vectors have different lengths."
	@assert k >= 1 "Overlap must be at least 1."
	@assert K >= 1  "Number of segments for approximating exp must be at least 1."
	@assert T <= length(q)  "Number of tests cannot exceed number of people in population."
	@assert all(isinteger.(u)) && all(u .>= 0)  "Utilities must be non-negative integers."
	@assert all(0 .<= q .<= 1)  "Probabilities must lie strictly between 0 and 1."

	# Compute some constants
	C = length(q)  # population size
	pset = collect(powerset(1:T, 1, k))  # array of subsets of tests of cardinality 1 to k
	evensets = [S for S in pset if iseven(length(S))]
	oddsets = [S for S in pset if isodd(length(S))]
	# println("pset: $(pset)")
	# println("evensets: $(evensets)")
	# println("oddsets: $(oddsets)")
	L, U = minimum(u), k*G*maximum(u)
	# println("L: $(L), U:$(U)")
	zvals = setdiff(L:U, 0)
	# println("zvals: $(zvals)")
	A, B = minimum(log.(u)) + k*G*minimum(log.(q)), log(k*G*maximum(u)) + maximum(log.(q))
	# println("A: $(A), B:$(B)")

	# Create model and set parameters
	m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV)))
	# "Presolve" => 0, "OutputFlag" => 0, "MIPGap" => 0.01, "TimeLimit" => 600
	!verbose && set_silent(m)
	
	@variable(m, x[1:T, 1:C], binary=true)  # test vectors
	@variable(m, con[pset, 1:C], binary=true)  # conjunctive variables
	@variable(m, dis[pset, 1:C], binary=true)  # disjunctive variables
	@variable(m, w[pset] >= 0)  # welfare for each set S in pset
	@variable(m, l[pset])  # log welfare for each set S in pset
	@variable(m, y[pset])
	@variable(m, z[pset])
	@variable(m, δ[oddsets, 1:K], binary=true)  # indicator variables for exp constraint and odd sets
	@variable(m, v[oddsets, 1:K])  # variables for exp constraint and odd sets
	@variable(m, γ[pset, union(zvals,0)], binary=true)  # indicator variables for log constraint

	@objective(m, Max, sum((-1)^(length(S)+1)*w[S] for S in pset))

	for S in pset  # ensure properties of conjunctive and disjunctive variables
		@constraint(m, [s in S, i in 1:C], con[S,i] <= x[s,i])  # Conjunctive
		@constraint(m, [i in 1:C], con[S,i] >= 1-length(S)+sum(x[s,i] for s in S))
		@constraint(m, [s in S, i in 1:C], dis[S,i] >= x[s,i])  # Disjunctive
		@constraint(m, [i in 1:C], dis[S,i] <= sum(x[s,i] for s in S))
	end

	# Deal with exponential constraints for even S
	a, b, c, ε = lower_bounds(A, B, 10*K)  # compute 10*K+1 linear lower bounds for exp on domain (A,B)
	error = length(evensets)*ε
	@constraint(m, [S in evensets, k in 1:10*K+1], w[S] >= a[k]*l[S]+b[k])
	# println("a: $(a), b: $(b), c:$(c), error:$(length(evensets)*ε)")
	# Deal with exponential constraints for odd S
	d, e, f, ε = upper_bounds(A, B, K)  # compute piecewise-linear function f on domain [A, B]
	error += length(oddsets)*ε
	# println("d: $(d), e: $(e), f:$(f), error: $(length(oddsets)*ε)")
	@constraint(m, c4b[S in oddsets], w[S] <= sum(d[k]*v[S,k] + e[k]*δ[S,k] for k in 1:K))
	@constraint(m, c4c[S in oddsets], sum(δ[S,k] for k in 1:K) <= 1)
	@constraint(m, c4d[S in oddsets], l[S] == sum(v[S,k] for k in 1:K))
	@constraint(m, c4e[S in oddsets, k=1:K], f[k]*δ[S,k] <= v[S,k])
	@constraint(m, c4f[S in oddsets, k=1:K], v[S,k] <= f[k+1]*δ[S,k])
	# Exponential constraints done
	@constraint(m, c4g[S in pset], l[S] == y[S] + sum(dis[S,i] * log(q[i]) for i in 1:C))
	# Start logarithmic constraints
	@constraint(m, c4h[S in pset], 1 == γ[S,0] + sum(γ[S,k] for k in zvals))
	@constraint(m, c4i[S in pset], z[S] == sum(k*γ[S,k] for k in zvals))
	@constraint(m, c4j[S in pset], y[S] == -1e6*γ[S,0] + sum(log(k)*γ[S,k] for k in zvals))
	# Logarithmic constraints done
	@constraint(m, c4k[S in pset], z[S] == sum(con[S,i] * u[i] for i in 1:C))  # utility sums
	@constraint(m, c4l[i=1:C], sum(x[t,i] for t in 1:T) <= k)  # k-overlapping tests
	@constraint(m, c4mn[t=1:T], 1 <= sum(x[t,i] for i in 1:C) <= G)  # pool size between 1 and G
	# Return model
	# return m, x, con, dis, w, l, y, z, γ, δ, v
	return m, x, error
end