@testset "bisect_interval" begin
    a, b = Arf(0), Arf(2)
    @test isequal(ArbExtras.bisect_interval(a, b), ((a, Arf(1)), (Arf(1), b)))
    a, b = Arf(0), Arf(Inf)
    @test isequal(ArbExtras.bisect_interval(a, b), ((a, Arf(Inf)), (Arf(Inf), b)))
    a, b, = Arf(0, prec = 80), Arf(2, prec = 90)
    @test precision(ArbExtras.bisect_interval(a, b)[1][2]) == 90

    a, b = Arf(1), Arf(2)
    @test ArbExtras.bisect_interval(a, b, log_midpoint = true)[1][2] ≈ exp(log(2) / 2)
    a, b = Arf(-2), Arf(-1)
    @test ArbExtras.bisect_interval(a, b, log_midpoint = true)[1][2] ≈ -exp(log(2) / 2)
    a, b = Arf(0), Arf(2)
    @test isequal(
        ArbExtras.bisect_interval(a, b, log_midpoint = true),
        ((a, Arf(1)), (Arf(1), b)),
    )
    a, b = Arf(-2), Arf(0)
    @test isequal(
        ArbExtras.bisect_interval(a, b, log_midpoint = true),
        ((a, Arf(-1)), (Arf(-1), b)),
    )
    a, b = Arf(-2), Arf(2)
    @test isequal(
        ArbExtras.bisect_interval(a, b, log_midpoint = true),
        ((a, Arf(0)), (Arf(0), b)),
    )

    a, b = Arb(0), Arb(2)
    @test isequal(ArbExtras.bisect_interval(a, b), ((a, Arb(1)), (Arb(1), b)))
    a, b = Arb(0), Arb(Inf)
    @test isequal(ArbExtras.bisect_interval(a, b), ((a, Arb(Inf)), (Arb(Inf), b)))
    a, b, = Arb(0, prec = 80), Arb(2, prec = 90)
    @test precision(ArbExtras.bisect_interval(a, b)[1][2]) == 90
end

@testset "bisect_interval_recursive" begin
    a, b = Arf(0), Arf(8)
    @test ArbExtras.bisect_interval_recursive(a, b, 0) == [(0, 8)]
    @test ArbExtras.bisect_interval_recursive(a, b, 1) == [(0, 4), (4, 8)]
    @test ArbExtras.bisect_interval_recursive(a, b, 2) == [(0, 2), (2, 4), (4, 6), (6, 8)]
    @test ArbExtras.bisect_interval_recursive(a, b, 3) ==
          [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8)]

    a, b = Arb(0), Arb(8)
    @test ArbExtras.bisect_interval_recursive(a, b, 0) == [(0, 8)]
    @test ArbExtras.bisect_interval_recursive(a, b, 1) == [(0, 4), (4, 8)]
    @test ArbExtras.bisect_interval_recursive(a, b, 2) == [(0, 2), (2, 4), (4, 6), (6, 8)]
    @test ArbExtras.bisect_interval_recursive(a, b, 3) ==
          [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8)]
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
