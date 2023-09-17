@testset "utilities" begin
    @testset "bisect_interval: $T" for T in (Arf, Arb)
        @test ArbExtras.bisect_interval(T(0), T(2)) == ((0, 1), (1, 2))
        @test ArbExtras.bisect_interval(T(0), T(Inf)) == ((0, Inf), (Inf, Inf))
        @test precision(
            ArbExtras.bisect_interval(T(0, prec = 80), T(2, prec = 90))[1][2],
        ) == 90

        @test ArbExtras.bisect_interval(T(1), T(4), log_midpoint = true) == ((1, 2), (2, 4))
        @test ArbExtras.bisect_interval(T(-4), T(-1), log_midpoint = true) ==
              ((-4, -2), (-2, -1))
        @test ArbExtras.bisect_interval(T(0), T(2), log_midpoint = true) == ((0, 1), (1, 2))
        @test ArbExtras.bisect_interval(T(-2), T(0), log_midpoint = true) ==
              ((-2, -1), (-1, 0))

        @test ArbExtras.bisect_interval(T(-2), T(2), log_midpoint = true) ==
              ((-2, 0), (0, 2))
    end

    @testset "bisect_interval_recursive: $T" for T in (Arf, Arb)
        a, b = T(0), T(8)
        @test ArbExtras.bisect_interval_recursive(a, b, 0) == [(0, 8)]
        @test ArbExtras.bisect_interval_recursive(a, b, 1) == [(0, 4), (4, 8)]
        @test ArbExtras.bisect_interval_recursive(a, b, 2) ==
              [(0, 2), (2, 4), (4, 6), (6, 8)]
        @test ArbExtras.bisect_interval_recursive(a, b, 3) ==
              [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8)]

        @test ArbExtras.bisect_interval_recursive(a, b, 0, log_midpoint = true) == [(0, 8)]
        @test ArbExtras.bisect_interval_recursive(a, b, 1, log_midpoint = true) ==
              [(0, 4), (4, 8)]
        @test all(
            map(
                ((x1, y1), (x2, y2)) -> x1 ≈ x2 && y1 ≈ y2,
                ArbExtras.bisect_interval_recursive(a, b, 2, log_midpoint = true),
                [(0, 2), (2, 4), (4, sqrt(32)), (sqrt(32), 8)],
            ),
        )
    end

    @testset "bisect_intervals: $T" for T in (Arf, Arb)
        intervals = [(T(0), T(2))]
        @test ArbExtras.bisect_intervals(intervals, [true]) ==
              ArbExtras.bisect_intervals(intervals, BitVector((true,))) ==
              [(0, 1), (1, 2)]
        @test ArbExtras.bisect_intervals(intervals, [false]) ==
              ArbExtras.bisect_intervals(intervals, BitVector((false,))) ==
              []

        intervals = [(T(0), T(2)), (T(4), T(6))]
        @test ArbExtras.bisect_intervals(intervals, [true, true]) ==
              ArbExtras.bisect_intervals(intervals, BitVector((true, true))) ==
              [(0, 1), (1, 2), (4, 5), (5, 6)]
        @test ArbExtras.bisect_intervals(intervals, [true, false]) ==
              ArbExtras.bisect_intervals(intervals, BitVector((true, false))) ==
              [(0, 1), (1, 2)]
        @test ArbExtras.bisect_intervals(intervals, [false, true]) ==
              ArbExtras.bisect_intervals(intervals, BitVector((false, true))) ==
              [(4, 5), (5, 6)]
        @test ArbExtras.bisect_intervals(intervals, [false, false]) ==
              ArbExtras.bisect_intervals(intervals, BitVector((false, false))) ==
              []

        intervals = [(T(1), T(4))]
        @test ArbExtras.bisect_intervals(intervals, [true], log_midpoint = true) ==
              ArbExtras.bisect_intervals(
                  intervals,
                  BitVector((true,)),
                  log_midpoint = true,
              ) ==
              [(1, 2), (2, 4)]
        @test ArbExtras.bisect_intervals(intervals, [false], log_midpoint = true) ==
              ArbExtras.bisect_intervals(
                  intervals,
                  BitVector((false,)),
                  log_midpoint = true,
              ) ==
              []
    end

    @testset "check_tolerance" begin
        @testset "x::Arb" begin
            @test ArbExtras.check_tolerance(Arb((1, 1.1)), atol = 0.2, rtol = 0.2)
            @test ArbExtras.check_tolerance(Arb((1, 1.1)), atol = 0.2, rtol = 0.05)
            @test ArbExtras.check_tolerance(Arb((1, 1.1)), atol = 0.05, rtol = 0.2)
            @test !ArbExtras.check_tolerance(Arb((1, 1.1)), atol = 0.05, rtol = 0.05)

            @test ArbExtras.check_tolerance(Arb((3, 3.1)), atol = 0.2)
            @test !ArbExtras.check_tolerance(Arb((3, 3.1)), atol = 0.05)
            @test ArbExtras.check_tolerance(Arb(1), atol = 0)

            @test ArbExtras.check_tolerance(Arb((10, 10.1)), rtol = 0.02)
            @test !ArbExtras.check_tolerance(Arb((10, 10.1)), rtol = 0.005)
            @test ArbExtras.check_tolerance(Arb(1), rtol = 0)
            @test ArbExtras.check_tolerance(Arb(0), rtol = 0)
            @test !ArbExtras.check_tolerance(Arb(π) - Arb(π), rtol = 1)

            @test !ArbExtras.check_tolerance(Arb(NaN))
            @test !ArbExtras.check_tolerance(Arb(Inf))
            @test ArbExtras.check_tolerance(Arb((-1000, 1000)))
        end

        @testset "x::AbstractVector" begin
            prop_check_tolerance_elementwise((x, atol, rtol)) =
                ArbExtras.check_tolerance(x; atol, rtol) ==
                all(xᵢ -> ArbExtras.check_tolerance(xᵢ; atol, rtol), x)

            sample = [
                Arb(0),
                Arb(1),
                Arb(π),
                Arb((1, 1.1)),
                Arb((3, 3.1)),
                Arb((10, 10.1)),
                Arb(NaN),
                Arb(Inf),
                Arb((-1000, 1000)),
            ]

            gen_tol = isample(range(0, 1, 10), PropCheck.noshrink)
            gen_vec = PropCheck.vector(isample(0:5), isample(sample, PropCheck.noshrink))

            gen = interleave(gen_vec, gen_tol, gen_tol)

            @test check(prop_check_tolerance_elementwise, gen, ntests = 1000)
        end
    end

    @testset "check_interval" begin
        @test ArbExtras.check_interval(Bool, Arf(1), Arf(2))
        @test ArbExtras.check_interval(Bool, Arf(-5), Arf(5))
        @test ArbExtras.check_interval(Bool, Arf(1), Arf(1))
        @test !ArbExtras.check_interval(Bool, Arf(2), Arf(1))
        @test !ArbExtras.check_interval(Bool, Arf(1), Arf(Inf))
        @test !ArbExtras.check_interval(Bool, Arf(-Inf), Arf(2))
        @test !ArbExtras.check_interval(Bool, Arf(-Inf), Arf(Inf))
        @test !ArbExtras.check_interval(Bool, Arf(Inf), Arf(2))
        @test !ArbExtras.check_interval(Bool, Arf(1), Arf(-Inf))
        @test !ArbExtras.check_interval(Bool, Arf(1), Arf(NaN))
        @test !ArbExtras.check_interval(Bool, Arf(NaN), Arf(2))

        @test isnothing(ArbExtras.check_interval(Arf(1), Arf(2)))
        @test isnothing(ArbExtras.check_interval(Arf(-5), Arf(5)))
        @test isnothing(ArbExtras.check_interval(Arf(1), Arf(1)))
        @test_throws ArgumentError ArbExtras.check_interval(Arf(2), Arf(1))
        @test_throws ArgumentError ArbExtras.check_interval(Arf(1), Arf(Inf))
        @test_throws ArgumentError ArbExtras.check_interval(Arf(-Inf), Arf(2))
        @test_throws ArgumentError ArbExtras.check_interval(Arf(-Inf), Arf(Inf))
        @test_throws ArgumentError ArbExtras.check_interval(Arf(Inf), Arf(2))
        @test_throws ArgumentError ArbExtras.check_interval(Arf(1), Arf(-Inf))
        @test_throws ArgumentError ArbExtras.check_interval(Arf(1), Arf(NaN))
        @test_throws ArgumentError ArbExtras.check_interval(Arf(NaN), Arf(2))
    end

    @testset "format_interval" begin
        @test ArbExtras.format_interval(Arf(1), Arf(2)) == "[1.0, 2.0]"
        @test ArbExtras.format_interval(Arf(1), Arf(Inf)) == "[1.0, Inf]"
        @test ArbExtras.format_interval(Arf(-Inf), Arf(2)) == "[-Inf, 2.0]"
        @test ArbExtras.format_interval(Arf(-Inf), Arf(Inf)) == "[-Inf, Inf]"

        @test ArbExtras.format_interval(Arf(3.1415), ubound(Arb(π))) == "[3.1415, 3.14159]"
        @test ArbExtras.format_interval(Arf(3.14159), ubound(Arb(π))) ==
              "[3.14159 +/- 2.66e-6]"
        @test ArbExtras.format_interval(lbound(Arb(π)), Arf(Inf)) == "[3.14159, Inf]"
        @test ArbExtras.format_interval(Arf(0), ubound(Arb(π))) == "[0.0, 3.14159]"

        @test ArbExtras.format_interval(Arf(1.23456789), Arf(1.234567891), digits = 10) ==
              "[1.23456789, 1.234567891]"
    end

    @testset "taylor_remainder" begin
        @test let x = Arb((-1, 1))
            Arblib.overlaps(
                ArbExtras.taylor_remainder(cos(ArbSeries((x, 1), degree = 4)), x),
                Arb(1 // factorial(4)),
            )
        end

        for f in (cos, exp, atan)
            for x in (Arb((1, 2)), Arb(NaN), Arb(5), Arb((-10, 10)))
                for degree = 1:10
                    p = f(ArbSeries((x, 1); degree))
                    @test isequal(
                        ArbExtras.taylor_remainder(p, x),
                        p[degree] * (x - midpoint(x))^degree,
                    )
                end
            end
        end
    end

    @testset "enclosure endpoints" begin
        x = Arb(3, prec = 16)
        r = Mag(lbound(Arb(π, prec = 16) - x))
        Arblib.set!(Arblib.radref(x), r)

        x = Arb(0, prec = 16)
        r = Mag(lbound(Arb(π, prec = 16)))
        Arblib.set!(Arblib.radref(x), r)

        @test -Arb(π, prec = 256) < x < Arb(π, prec = 256)

        @test Arb(π, prec = 256) < ubound(Arb, x)
        @test lbound(Arb, x) < -Arb(π, prec = 256)

        @test Arblib.overlaps(ArbExtras.enclosure_ubound(x), Arb(π, prec = 256))
        @test Arblib.overlaps(ArbExtras.enclosure_lbound(x), -Arb(π, prec = 256))

        @test isequal(
            (ArbExtras.enclosure_lbound(x), ArbExtras.enclosure_ubound(x)),
            ArbExtras.enclosure_getinterval(x),
        )

        for x in Arb.([1, 1 // 3, ℯ, π, (0, 1)])
            @test Arblib.overlaps(ubound(Arb, x), ArbExtras.enclosure_ubound(x))
            @test Arblib.overlaps(lbound(Arb, x), ArbExtras.enclosure_lbound(x))

            @test isequal(
                (ArbExtras.enclosure_lbound(x), ArbExtras.enclosure_ubound(x)),
                ArbExtras.enclosure_getinterval(x),
            )
        end
    end

    @testset "derivative_function" begin
        args = (
            Arb(0.5),
            Arb(1),
            Arb(1.5),
            Acb(0.5, 0.6),
            Acb(1, 1.1),
            Acb(1.5, 1.6),
            ArbSeries((0.5, 1)),
            ArbSeries((1, 2.5)),
            ArbSeries((1.5, -1)),
            ArbSeries((1, 2, 3)),
            ArbSeries((1, 2, 3, 4)),
            ArbSeries((-1, 2, -3, 4)),
            ArbSeries((0.5, 1), degree = 5),
            AcbSeries((0.5 + 0.6im, 1)),
            AcbSeries((1 + 1.1im, 2.5 + 2.6im)),
            AcbSeries((1.5 + 1.6im, -1)),
            AcbSeries((1 + 1.1im, 2 + 2.2im, 3 + 3.3im)),
            AcbSeries((1 + 1.1im, 2 + 2.2im, 3 + 3.3im, 4 + 4.4im)),
            AcbSeries((-1 + 1.1im, 2 - 2.2im, -3 + 3.3im, 4 - 4.4im)),
            AcbSeries((0.5 + 0.6im, 1), degree = 5),
        )

        for (f1, df1, d2f1, d3f1) in (
            (sin, cos, x -> -sin(x), x -> -cos(x)),
            (
                tan,
                x -> sec(x)^2,
                x -> 2tan(x) * sec(x)^2,
                x -> 2sec(x)^2 * (2tan(x)^2 + sec(x)^2),
            ),
            (
                atan,
                x -> inv(1 + x^2),
                x -> -2x / (1 + x^2)^2,
                x -> (6x^2 - 2) / (1 + x^2)^3,
            ),
        )

            f2 = ArbExtras.derivative_function(f1, 0)
            df2 = ArbExtras.derivative_function(f1, 1)
            d2f2 = ArbExtras.derivative_function(f1, 2)
            d3f2 = ArbExtras.derivative_function(f1, 3)

            for x in args
                @test Arblib.overlaps(f1(x), f2(x))
                @test Arblib.overlaps(df1(x), df2(x))
                @test Arblib.overlaps(d2f1(x), d2f2(x))
                @test Arblib.overlaps(d3f1(x), d3f2(x))
            end
        end
    end
end
