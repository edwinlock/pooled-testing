@testset "solvers (Gurobi + MOSEK)" begin
    if !SOLVERS_AVAILABLE
        @info "Skipping solver tests — Gurobi/MOSEK licence not available."
    else
        # A small fixed population for deterministic checks.
        pop = Population{Int}(1 => (0.9, 10), 2 => (0.85, 8), 3 => (0.8, 6),
                              4 => (0.7, 5), 5 => (0.6, 4))

        @testset "greedy" begin
            w, pools, err = greedy(pop; T=2, G=3)
            @test w > 0
            @test err == 0.0                              # greedy has no approximation error
            @test length(pools) ≤ 2                       # at most T pools
            @test welfare(pools, pop) ≈ w                 # reported welfare matches the pools
        end

        @testset "exact (disjoint, k=1)" begin
            w, pools, err = exact(pop; k=1, T=2, G=3)
            @test w > 0
            @test welfare(pools, pop) ≈ w
            # Pools are disjoint (non-overlapping): no individual appears twice.
            members = reduce(vcat, pools; init=Int[])
            @test length(members) == length(unique(members))
        end

        @testset "approximate honours the (welfare, pools, error) contract" begin
            w, pools, err = approximate(pop; T=2, G=3, K=15)
            @test w > 0
            @test err ≥ 0                                 # additive guarantee, non-negative
            @test welfare(pools, pop) ≈ w atol=1e-6
        end

        @testset "approximate ≈ exact within the guarantee" begin
            wa, _, err = approximate(pop; T=2, G=3, K=25)
            we, _, _   = exact(pop; k=1, T=2, G=3)
            # The exact non-overlapping optimum is within the MILP's additive guarantee.
            @test we ≤ wa + err + 1e-6
            @test wa ≤ we + err + 1e-6
        end

        @testset "empty population returns the 3-tuple contract" begin
            empty_pop = Population{Int}(1 => (0.0, 0))     # removed by remove_zeros!
            @test approximate(empty_pop; T=1, G=2) == (0.0, [], 0.0)
            @test exact(empty_pop; k=1, T=1, G=2) == (0.0, [], 0.0)
        end
    end
end
