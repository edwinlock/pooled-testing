using MosekTools, MathOptInterface
const MOI = MathOptInterface


"""
MICOP model for a single test that can be solved exactly with MOSEK.
	
Note that the objective value of this model is the log of the overall welfare!
"""
function single_mosek(q, u; G, verbose=false)
	# Verify that input is consistent
	@assert length(u) == length(q) "Input vectors have different lengths."
	@assert all(u .>= 0) "Utilities must be non-negative."
	@assert all(0 .<= q .<= 1) "Probabilities must lie between 0 and 1."
	
	# Compute population size
	n = length(u)

	# Create model and set parameters
	m = Model(Mosek.Optimizer)
	# Force MOSEK to a single internal thread. The outer experiment loop already
	# parallelises across solves; letting MOSEK spin up its own (TBB) thread pool
	# from within a Julia worker thread triggers a TBB init crash on many-core
	# machines. One thread per solve avoids that and is best for these small MICPs.
	set_optimizer_attribute(m, "MSK_IPAR_NUM_THREADS", 1)
	!verbose && set_silent(m)

	# Define variables
	@variable(m, x[[1],1:n], binary=true)
	@variable(m, l)
	@variable(m, y)
	@variable(m, z)

	# Add objective
	@objective(m, Max, l)

	# Add constraints
	# z = u ̇x and y <= log(u ̇x)
	@constraint(m, z == sum(x[1,i]*u[i] for i in 1:n))
	@constraint(m, [y, 1, z] in MOI.ExponentialCone())
	# Welfare constraint: l <= log(u ̇x) + x ̇log.(q)
	@constraint(m, l == y + sum(x[1,i]*log(q[i]) for i in 1:n))
	# Pooled testing size constraint
	@constraint(m, sum(x[1,i] for i in 1:n) <= G)

	# Return model
	return m, x
end
