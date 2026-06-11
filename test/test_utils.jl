@testset "population utilities" begin

    @testset "welfare (inclusion-exclusion)" begin
        pop = Population{Int}(1 => (0.9, 10), 2 => (0.8, 5), 3 => (0.5, 8))
        # Single pool {1,2}: (u1+u2)·q1·q2 = 15·0.72 = 10.8
        @test welfare([[1, 2]], pop) ≈ 10.8
        # Singleton {1}: u1·q1 = 10·0.9
        @test welfare([[1]], pop) ≈ 9.0
        # Two disjoint pools add up (no shared individuals): 10.8 + 8·0.5
        @test welfare([[1, 2], [3]], pop) ≈ 10.8 + 4.0
        # Empty allocation gives zero welfare.
        @test welfare([], pop) == 0
    end

    @testset "generate_instance" begin
        pop = generate_instance(7, 0:0.1:1, 1:10)
        @test length(pop) == 7
        @test all(0 <= v[1] <= 1 for v in values(pop))   # probabilities in range
        @test all(1 <= v[2] <= 10 for v in values(pop))  # utilities in range
        @test all(isinteger(v[2]) for v in values(pop))  # utilities rounded to Int
    end

    @testset "scale_utilities" begin
        # Below the cap: utilities unchanged.
        pop = Population{Int}(1 => (0.9, 10), 2 => (0.5, 30))
        @test scale_utilities(pop, 50) == pop
        # Above the cap: scaled into [1, upper], max maps to upper.
        big = Population{Int}(1 => (0.9, 1000), 2 => (0.8, 500))
        s = scale_utilities(big, 50)
        @test all(1 <= v[2] <= 50 for v in values(s))
        @test s[1][2] == 50
        # Tiny utilities floor at 1 rather than rounding to 0 (no one dropped).
        tiny = Population{Int}(1 => (0.9, 1000), 2 => (0.8, 1))
        @test scale_utilities(tiny, 50)[2][2] == 1
        # Empty population is returned untouched.
        @test isempty(scale_utilities(Population{Int}(), 50))
    end

    @testset "remove_zeros!" begin
        pop = Population{Int}(1 => (0.9, 10), 2 => (0.0, 5), 3 => (0.7, 0))
        remove_zeros!(pop)
        @test collect(keys(pop)) == [1]   # drops zero-prob (2) and zero-util (3)
    end

    @testset "pop2vec round-trips" begin
        pop = Population{Int}(1 => (0.9, 10), 2 => (0.8, 5))
        q, u, keys = pop2vec(pop)
        @test length(q) == length(u) == length(keys) == 2
        for (i, k) in enumerate(keys)            # vectors align with keylist
            @test (q[i], u[i]) == pop[k]
        end
    end

    @testset "clustering" begin
        # Two individuals share a (q,u); a third differs → two clusters, sizes 2 and 1.
        pop = Population{Int}(1 => (0.9, 10), 2 => (0.9, 10), 3 => (0.8, 5))
        clusters = pop2clusters(pop)
        @test sort(length.(values(clusters))) == [1, 2]
        q, u, n, keylist = cluster2vec(clusters)
        @test sort(n) == [1, 2]
        @test sum(n) == 3                        # cluster sizes cover everyone
        @test Set(keylist) == Set([(0.9, 10), (0.8, 5)])
    end

    @testset "uncluster maps indices back to people" begin
        pop = Population{Int}(1 => (0.9, 10), 2 => (0.9, 10), 3 => (0.8, 5))
        clusters = pop2clusters(pop)
        _, _, _, keylist = cluster2vec(clusters)
        # Pool referencing cluster index 1 twice draws two distinct members of it.
        i = findfirst(==( (0.9, 10) ), keylist)
        unclustered = uncluster([[i, i]], deepcopy(clusters), keylist)
        @test length(unclustered[1]) == 2
        @test Set(unclustered[1]) == Set([1, 2])  # the two members of that cluster
    end

end
