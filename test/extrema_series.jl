@testset "extrema_series" begin
    @testset "test problems" begin
        demo_problems = demo_problem.(eachindex(demo_problems_definitions), true)
        for (f, (a, b), min) in demo_problems
            min1, _, y1 = ArbExtras.extrema_series(f, a, b)
            _, min2, y2 = .-ArbExtras.extrema_series(x -> -f(x), a, b)
            min3, y3 = ArbExtras.minimum_series(f, a, b)
            min4, y4 = .-ArbExtras.maximum_series(x -> -f(x), a, b)

            @test isequal(min1, min2)
            @test isequal(min2, min3)
            @test isequal(min3, min4)

            @test isequal(y1, y2)
            @test isequal(y3, y3)
            @test isequal(y3, y4)

            @test Arblib.overlaps(min1, min)

            @test Arblib.overlaps(y1, f(Arb(a + b) / 2))
        end
    end

    @testset "special cases" begin
        # Thin interval
        @test ArbExtras.extrema_series(x -> -2x, Arf(2), Arf(2)) == (-4, -4, -4)
        @test ArbExtras.minimum_series(x -> -2x, Arf(2), Arf(2)) == (-4, -4)
        @test ArbExtras.maximum_series(x -> -2x, Arf(2), Arf(2)) == (-4, -4)
        @test ArbExtras.extrema_series(x -> -2x, Arf(2), Arf(2), abs_value = true) ==
              (4, 4, 4)
        @test ArbExtras.minimum_series(x -> -2x, Arf(2), Arf(2), abs_value = true) == (4, 4)
        @test ArbExtras.maximum_series(x -> -2x, Arf(2), Arf(2), abs_value = true) == (4, 4)

        # Non-ordered endpoints
        @test_throws ArgumentError ArbExtras.extrema_series(identity, Arf(2), Arf(1))
        @test_throws ArgumentError ArbExtras.minimum_series(identity, Arf(2), Arf(1))
        @test_throws ArgumentError ArbExtras.maximum_series(identity, Arf(2), Arf(1))

        # Infinite endpoints
        @test_throws ArgumentError ArbExtras.extrema_series(identity, Arf(2), Arf(Inf))
        @test_throws ArgumentError ArbExtras.minimum_series(identity, Arf(2), Arf(Inf))
        @test_throws ArgumentError ArbExtras.maximum_series(identity, Arf(2), Arf(Inf))

        # Monotonic on interval
        @test isequal(ArbExtras.extrema_series(identity, Arf(-2), Arf(4)), (-2, 4, NaN))
        @test isequal(ArbExtras.minimum_series(identity, Arf(-2), Arf(4)), (-2, NaN))
        @test isequal(ArbExtras.maximum_series(identity, Arf(-2), Arf(4)), (4, NaN))

        @test isequal(
            ArbExtras.extrema_series(identity, Arf(-2), Arf(4), abs_value = true),
            (0, 4, NaN),
        )
        @test isequal(
            ArbExtras.extrema_series(identity, Arf(-3), Arf(-1), abs_value = true),
            (1, 3, NaN),
        )
        @test isequal(
            ArbExtras.minimum_series(identity, Arf(-2), Arf(4), abs_value = true),
            (0, NaN),
        )
        @test isequal(
            ArbExtras.minimum_series(identity, Arf(-3), Arf(-1), abs_value = true),
            (1, NaN),
        )
        @test isequal(
            ArbExtras.maximum_series(identity, Arf(-6), Arf(4), abs_value = true),
            (6, NaN),
        )

        # Check that it throws an error if the computations with the
        # endpoints can not be done exactly
        @test_throws AssertionError ArbExtras.extrema_series(sin, Arf(1), Arf(1e300))
        @test_throws AssertionError ArbExtras.minimum_series(sin, Arf(1), Arf(1e300))
        @test_throws AssertionError ArbExtras.maximum_series(sin, Arf(1), Arf(1e300))

        # Non-finite restterm
        @test all(isnan, ArbExtras.extrema_series(x -> inv(x - sin(x)), Arf(0.1), Arf(1)))
        @test all(isnan, ArbExtras.minimum_series(x -> inv(x - sin(x)), Arf(0.1), Arf(1)))
        @test all(isnan, ArbExtras.maximum_series(x -> inv(x - sin(x)), Arf(0.1), Arf(1)))
    end
end
