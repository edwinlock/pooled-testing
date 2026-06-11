# Test suite for the pooled-testing code.
#
#   julia --project=. test/runtests.jl
#
# Solver-free tests (population utilities, the piecewise-linear approximation
# math, the SQLite store, analysis reshaping) always run. Tests that need Gurobi
# and MOSEK licences run only when a licence is available; otherwise they are
# skipped with a notice, so the suite passes anywhere (e.g. CI).

include("test_setup.jl")  # loads the project code once; sets SOLVERS_AVAILABLE

@testset "pooled-testing" begin
    include("test_utils.jl")            # population helpers, welfare, scaling, clustering
    include("test_approximation.jl")    # piecewise-linear approximation + guarantee
    include("test_datastore.jl")        # SQLite store: hashing, keys, idempotency, reuse
    include("test_analyse.jl")          # reshaping the store into summary frames
    include("test_solvers.jl")          # greedy / exact / approximate (licence-gated)
end
