using JuMP, Ipopt, SCIP, Combinatorics # Alpine, Gurobi, Juniper, SCIP

# include("optimizers.jl")

# # Choose underlying solvers for Alpine
# nlp_solver = get_ipopt() # local continuous solver
# mip_solver = get_gurobi() # convex mip solver
# minlp_solver = get_juniper(mip_solver, nlp_solver) # local mixed-intger solver

# const alpine = JuMP.optimizer_with_attributes(
#     Alpine.Optimizer,
#     "nlp_solver" => nlp_solver,
#     "mip_solver" => mip_solver,
#     "minlp_solver" => minlp_solver,
#     "presolve_bt" => true,
#     "apply_partitioning" => true,
#     "partition_scaling_factor" => 10,
# )

function create_problem(n)
    m = Model(SCIP.Optimizer)

    # Compute all possible allocations using k ∈ {1,2,3} tests
    allocations1 = all_allocations(n, 1)
    allocations2 = all_allocations(n, 2)
    allocations3 = all_allocations(n, 3)
    alloc = [allocations1, allocations2, allocations3]
    all = union(allocations1, allocations2, allocations3)

    # Create variables
    @variable(m, u[1:n] >= 0)  # utilities
    @variable(m, 0 <= q[1:n] <= 1)  # probabilities
    @variable(m, w[1:3])  # maximum welfare achieveable for each budget
    @variable(m, welfare[all])  # welfare achieved by each possible allocation
    @variable(m, indicators[B in [1,3], alloc[B]], Bin)  # indicator variables to model the max constraints
    # NB: B stands for budget

    # Constraints
    @NLconstraint(  # welfare[a] should be the expected welfare of allocation `a` wrt u and q.
        m,
        [a in all],
        welfare[a] == sum(prod(q[i] for i in pool) * sum(u[i] for i in pool) for pool in a)
    )
    # Model constraint w[B] == max(welfares[a] for a in alloc[B]) for all B ∈ {1,2,3} with indicator variables:
    # First ensure w[B] >= w[a] for all a in alloc[B]
    @constraint(m, [B in 1:3, a in alloc[B]], w[B] >= welfare[a])
    # Now ensure that $w[B] == w[a] for some a, for B=1 and B=3 (B=2 is not necessary!)
    @constraint(m, [B in [1,3]], sum(indicators[B,a] for a in alloc[B]) == 1)
    @NLconstraint(m, [B in [1,3]], w[B] == sum(indicators[B, a] * welfare[a] for a in alloc[B]))

    # Objective
    @objective(m, Max, w[1] - 2*w[2] + w[3])
    
    return m

end

m = create_problem(3)
optimize!(m)


## CODE GRAVEYARD

# # Register custom max function for each allocation (technical hack)  
# f(y...) = max()
# for B in 1:3
#     JuMP.register(m, Symbol("max$(B)"), length(alloc[B]), f; autodiff=true)
# end

# # Add maximum constraint
# alloc_welfares = []
# for B in 1:3
#     push!(alloc_welfares, [welfare[a] for a in alloc[B]])
# end
# @NLconstraint(m, [B in 1:3], w[B] == max1(alloc_welfares[B]...))