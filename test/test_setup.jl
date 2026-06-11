# Shared test scaffolding: load the package and detect whether the Gurobi/MOSEK
# solvers are usable, so solver tests can be skipped cleanly when they are not.

using Test, PooledTesting, DataFrames
using PooledTesting: Δ, partition, upper_bounds, exp_domain, milp_guarantee,
    remove_zeros!, pop2vec, pop2clusters, cluster2vec, uncluster,
    GUROBI_MIPGAP, NON_SOLVE_ARGS, conic

"""
True if the Gurobi/MOSEK solvers actually run (licences present). Probed once by
attempting a trivial `approximate` solve; solver tests skip cleanly when absent.
"""
const SOLVERS_AVAILABLE = try
    probe = Population{Int}(1 => (0.9, 5), 2 => (0.8, 3))
    approximate(probe; T=1, G=2, K=5)   # Gurobi
    greedy(probe; T=1, G=2)             # MOSEK (via conic)
    true
catch err
    @warn "Gurobi/MOSEK not available — solver tests will be skipped." exception=err
    false
end
