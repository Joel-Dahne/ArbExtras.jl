@testset "bisect_interval: $T" for T in (Arf, Arb)
    @test ArbExtras.bisect_interval(T(0), T(2)) == ((0, 1), (1, 2))
    @test ArbExtras.bisect_interval(T(0), T(Inf)) == ((0, Inf), (Inf, Inf))
    @test precision(ArbExtras.bisect_interval(T(0, prec = 80), T(2, prec = 90))[1][2]) == 90

    @test ArbExtras.bisect_interval(T(1), T(4), log_midpoint = true) == ((1, 2), (2, 4))
    @test ArbExtras.bisect_interval(T(-4), T(-1), log_midpoint = true) ==
          ((-4, -2), (-2, -1))
    @test ArbExtras.bisect_interval(T(0), T(2), log_midpoint = true) == ((0, 1), (1, 2))
    @test ArbExtras.bisect_interval(T(-2), T(0), log_midpoint = true) == ((-2, -1), (-1, 0))

    @test ArbExtras.bisect_interval(T(-2), T(2), log_midpoint = true) == ((-2, 0), (0, 2))
end

@testset "bisect_interval_recursive: $T" for T in (Arf, Arb)
    a, b = T(0), T(8)
    @test ArbExtras.bisect_interval_recursive(a, b, 0) == [(0, 8)]
    @test ArbExtras.bisect_interval_recursive(a, b, 1) == [(0, 4), (4, 8)]
    @test ArbExtras.bisect_interval_recursive(a, b, 2) == [(0, 2), (2, 4), (4, 6), (6, 8)]
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
          ArbExtras.bisect_intervals(intervals, BitVector((true,)), log_midpoint = true) ==
          [(1, 2), (2, 4)]
    @test ArbExtras.bisect_intervals(intervals, [false], log_midpoint = true) ==
          ArbExtras.bisect_intervals(intervals, BitVector((false,)), log_midpoint = true) ==
          []
end

@testset "check_tolerance" begin
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
    @test ArbExtras.format_interval(Arf(3.14159), ubound(Arb(π))) == "[3.14159 +/- 2.66e-6]"
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
