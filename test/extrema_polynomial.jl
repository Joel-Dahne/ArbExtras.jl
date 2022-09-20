@testset "extrema_polynomial" begin
    # Define a number of problems for which we expect the method to
    # work without issues. The problems are defined as (p, (a, b), x,
    # y). Where the minimum is attained at `x` and the maximum at `y`.
    problems = [
        # Min and max at interior points
        (
            Arblib.integral(Arblib.fromroots(ArbPoly, Arb[0.25, 1.75])),
            (Arf(0), Arf(2)),
            Arb(1.75),
            Arb(0.25),
        ),

        # Min at interior point, max at endpoint
        (ArbPoly([0, 0, 1]), (Arf(-1), Arf(2)), Arb(0), Arb(2)),
        (
            Arblib.integral(Arblib.fromroots(ArbPoly, [-1, 0, 1])),
            (Arf(-1.5), Arf(2)),
            Arb(1),
            Arb(2),
        ),
        (
            Arblib.integral(Arblib.fromroots(ArbPoly, [0, 1, 4])),
            (Arf(-1), Arf(6)),
            Arb(4),
            Arb(6),
        ),

        # Min at endpoint, max at interior point
        (ArbPoly([0, 0, -1]), (Arf(-1), Arf(2)), Arb(2), Arb(0)),
        (
            -Arblib.integral(Arblib.fromroots(ArbPoly, [-1, 0, 1])),
            (Arf(-1.5), Arf(2)),
            Arb(2),
            Arb(1),
        ),
        (
            -Arblib.integral(Arblib.fromroots(ArbPoly, [0, 1, 4])),
            (Arf(-1), Arf(6)),
            Arb(6),
            Arb(4),
        ),

        # Min and max at endpoints
        (ArbPoly([0, 1]), (Arf(-1), Arf(1)), Arb(-1), Arb(1)),
        (ArbPoly([1, -1, 3, -1]), (Arf(-3), Arf(3)), Arb(3), Arb(-3)),

        # Derivative has double root so cant isolate roots
        (
            Arblib.integral(Arblib.fromroots(ArbPoly, [0.25, 0.25, 1.75, 1.75])),
            (Arf(0), Arf(2)),
            Arb(0),
            Arb(2),
        ),

        # Special Cases

        # Thin interval
        (ArbPoly([1, 1, 1]), (Arf(1), Arf(1)), Arb(1), Arb(1)),

        # Constant polynomial
        (ArbPoly(5), (Arf(0), Arf(1)), Arb(0.5), Arb(0.5)),

        # Linear polynomial
        (ArbPoly([5, 1]), (Arf(0), Arf(1)), Arb(0), Arb(1)),
        (ArbPoly([5, -1]), (Arf(0), Arf(1)), Arb(1), Arb(0)),
    ]

    i = 0
    for (p, (a, b), x, y) in problems
        min, max = p(x), p(y)
        min1, max1 = ArbExtras.extrema_polynomial(p, a, b)
        min2 = ArbExtras.minimum_polynomial(p, a, b)
        max2 = ArbExtras.maximum_polynomial(p, a, b)

        @test isequal(min1, min2)
        @test isequal(max1, max2)

        @test Arblib.overlaps(min, min1)
        @test Arblib.overlaps(max, max1)

        @test radius(min1) < sqrt(eps(Arb))
        @test radius(max1) < sqrt(eps(Arb))
    end

    # Polynomial with NaN
    @test all(isnan, ArbExtras.extrema_polynomial(ArbPoly(NaN), Arf(0), Arf(1)))
    @test isnan(ArbExtras.minimum_polynomial(ArbPoly(NaN), Arf(0), Arf(1)))
    @test isnan(ArbExtras.maximum_polynomial(ArbPoly(NaN), Arf(0), Arf(1)))

    # TODO: Isolate roots fails
    #(Arblib.integral(Arblib.fromroots(ArbPoly, [0, 1, 2, 3, 4])), (Arf(-1), Arf(6)), Arb(0), Arb(6)),

    @testset "abs_value = true" begin
        # Define a number of problems for which we expect the method to
        # work without issues. The problems are defined as (p, (a, b), x,
        # y). Where the minimum is attained at `x` and the maximum at `y`.
        problems = [
            # Min and max at interior points
            (
                Arblib.integral(Arblib.fromroots(ArbPoly, [-1, 1])) + 2,
                (Arf(-1.5), Arf(1.5)),
                Arf(1),
                Arf(-1),
            ),
            (ArbPoly([-1, 0, 1]), (Arf(-1.1), Arf(1.2)), Arf(1), Arf(0)),

            # Min is non-zero at interior point, max at endpoint
            (ArbPoly([1, 0, 1]), (Arf(-1), Arf(2)), Arf(0), Arf(2)),
            (
                Arblib.integral(Arblib.fromroots(ArbPoly, [-0.5, 0, 0.5])) - 9 // 1024,
                (Arf(-0.9), Arf(1)),
                Arf(0.75),
                Arf(1),
            ),

            # Min is zero at interior point, max at endpoint
            (ArbPoly([-1, 0, 1]), (Arf(0), Arf(2)), Arf(1), Arf(2)),
            (ArbPoly([-1, 0, 1]), (Arf(-2), Arf(2)), Arf(1), Arf(2)),

            # Min at endpoint, max at interior point
            (ArbPoly([5, 0, -1]), (Arf(-1), Arf(2)), Arf(2), Arf(0)),

            # Min and max at endpoints
            (
                Arblib.integral(Arblib.fromroots(ArbPoly, [-1, 1])) + 2,
                (Arf(-2.2), Arf(2.2)),
                Arf(-2.2),
                Arf(2.2),
            ),

            # Derivative has double root so cant isolate roots
            (ArbPoly([5, 0, 0, 1]), (Arf(-1), Arf(2)), Arb(-1), Arb(2)),

            # Special Cases

            # Thin interval
            (ArbPoly([1, 1, 1]), (Arf(1), Arf(1)), Arb(1), Arb(1)),

            # Constant polynomial
            (ArbPoly(5), (Arf(0), Arf(1)), Arb(0.5), Arb(0.5)),

            # Linear polynomial
            (ArbPoly([0, 1]), (Arf(-2), Arf(1)), Arb(0), Arb(-2)),
            (ArbPoly([0, 1]), (Arf(-1), Arf(2)), Arb(0), Arb(2)),
            (ArbPoly([5, -1]), (Arf(0), Arf(1)), Arb(1), Arb(0)),
            (ArbPoly([-5, -1]), (Arf(0), Arf(1)), Arb(0), Arb(1)),

            # Zero at endpoint
            (ArbPoly([0, 1, 1]), (Arf(0), Arf(2)), Arb(0), Arb(2)),
        ]


        i = 0
        for (p, (a, b), x, y) in problems
            min, max = abs(p(x)), abs(p(y))
            # Run the problems with both p and -p, they should give
            # identical results
            min11, max11 = ArbExtras.extrema_polynomial(p, a, b, abs_value = true)
            min12, max12 = ArbExtras.extrema_polynomial(-p, a, b, abs_value = true)
            min21 = ArbExtras.minimum_polynomial(p, a, b, abs_value = true)
            min22 = ArbExtras.minimum_polynomial(-p, a, b, abs_value = true)
            max21 = ArbExtras.maximum_polynomial(p, a, b, abs_value = true)
            max22 = ArbExtras.maximum_polynomial(-p, a, b, abs_value = true)

            @test isequal(min11, min12)
            @test isequal(min11, min21)
            @test isequal(min11, min22)
            @test isequal(max11, max12)
            @test isequal(max11, max21)
            @test isequal(max11, max22)

            @test Arblib.overlaps(min, min11)
            @test Arblib.overlaps(max, max11)

            @test radius(min11) < sqrt(eps(Arb))
            @test radius(max11) < sqrt(eps(Arb))
        end

        # Polynomial with NaN
        @test all(
            isnan,
            ArbExtras.extrema_polynomial(ArbPoly(NaN), Arf(0), Arf(1), abs_value = true),
        )
        @test isnan(
            ArbExtras.minimum_polynomial(ArbPoly(NaN), Arf(0), Arf(1), abs_value = true),
        )
        @test isnan(
            ArbExtras.maximum_polynomial(ArbPoly(NaN), Arf(0), Arf(1), abs_value = true),
        )
    end
end
