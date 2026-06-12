@testset "piecewise-linear approximation" begin

    @testset "Δ: max gap between a chord of exp and exp itself" begin
        @test Δ(1.0, 1.0) == 0.0          # degenerate interval: no gap
        @test Δ(2.0, 1.0) == 0.0          # r ≤ l: defined as 0
        @test Δ(0.0, 1.0) > 0             # genuine chord lies above exp ⇒ positive gap
        # Wider interval ⇒ larger maximum gap.
        @test Δ(0.0, 2.0) > Δ(0.0, 1.0)
    end

    @testset "upper_bounds: K equal-error segments spanning [A,B]" begin
        A, B, K = 0.0, 5.0, 8
        a, b, c, ε = upper_bounds(A, B, K)
        @test length(c) == K + 1           # K+1 breakpoints
        @test c[1] ≈ A && c[end] ≈ B       # partition spans exactly [A,B]
        @test issorted(c)                  # breakpoints increasing
        @test length(a) == length(b) == K  # one (slope, intercept) per segment
        @test ε > 0                        # non-degenerate ⇒ positive error
        # Each segment over-approximates exp on its endpoints (f ≥ exp at breakpoints).
        for k in 1:K
            @test a[k]*c[k]   + b[k] ≥ exp(c[k])   - 1e-6
            @test a[k]*c[k+1] + b[k] ≥ exp(c[k+1]) - 1e-6
        end
    end

    @testset "error decreases as K grows (≈ 1/K²)" begin
        A, B = 0.0, 5.0
        err(K) = upper_bounds(A, B, K)[4]
        @test err(20) < err(10) < err(5)
        # Doubling K roughly quarters the error (piecewise-linear of a convex fn).
        @test err(10) / err(20) > 2        # at least a clear improvement, well above 1
    end

    @testset "exp_domain" begin
        q, u, G = [0.9, 0.8], [10, 5], 5
        A, B = exp_domain(q, u, G)
        @test A == minimum(log.(u)) + G*minimum(log.(q))
        @test B == log(G*maximum(u)) + maximum(log.(q))
        @test A < B
    end

    @testset "milp_guarantee = T·ε, scaling with T and K" begin
        q, u = [0.9, 0.8, 0.7], [10, 5, 8]
        g(T, K) = milp_guarantee(q, u; T=T, G=5, K=K)
        @test g(2, 15) > 0
        @test g(4, 15) ≈ 2 * g(2, 15)      # linear in budget T
        @test g(2, 20) < g(2, 10)          # tighter approximation with more segments
    end

    @testset "accuracy_params splits a total additive error budget" begin
        q, u = [0.9, 0.8, 0.7], [10, 5, 8]
        T, G, target, split = 2, 5, 0.25, 0.7
        K, gapabs = accuracy_params(q, u; T=T, G=G, target=target, split=split)
        formulation_error = milp_guarantee(q, u; T=T, G=G, K=K)

        @test K ≥ 1
        @test formulation_error ≤ split * target + 1e-9
        @test gapabs ≈ target - formulation_error atol=1e-9
        @test formulation_error + gapabs ≈ target atol=1e-9

        K2, gapabs2 = accuracy_params(q, u; T=T, G=G, target=target / 2, split=split)
        @test K2 ≥ K
        @test gapabs2 ≥ 0
        @test_throws AssertionError accuracy_params(q, u; T=T, G=G, target=0.0)
        @test_throws AssertionError accuracy_params(q, u; T=T, G=G, target=target, split=1.0)
    end

end
