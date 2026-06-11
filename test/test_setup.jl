# Shared test scaffolding: load the project code once, and detect whether the
# Gurobi/MOSEK solvers are usable so solver tests can be skipped when not.

using Test

# run.jl/analyse.jl define these before including optimisation.jl; tests do the same.
isdefined(Main, :GUROBI_MIPGAP)  || (const GUROBI_MIPGAP  = 1e-4)
isdefined(Main, :GUROBI_THREADS) || (const GUROBI_THREADS = 8)

const REPO = dirname(@__DIR__)
isdefined(Main, :welfare)     || include(joinpath(REPO, "optimisation.jl"))
isdefined(Main, :population_hash) || include(joinpath(REPO, "datastore.jl"))

"""
True if a Gurobi solve can actually run (licence present and working). Probed once
by attempting a trivial solve; solver tests use this to skip cleanly when absent.
"""
const SOLVERS_AVAILABLE = try
    using JuMP, Gurobi
    m = Model(() -> Gurobi.Optimizer(grb_env()))
    set_silent(m)
    @variable(m, 0 <= x <= 1)
    @objective(m, Max, x)
    optimize!(m)
    termination_status(m) == MOI.OPTIMAL
catch err
    @warn "Gurobi/MOSEK not available — solver tests will be skipped." exception=err
    false
end
