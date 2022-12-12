@testset "_check_root_interval" begin
    @testset "function" begin
        # Unique zero in the interval
        @test ArbExtras._check_root_interval(sin, Arf(-1), Arf(1)) == (true, true)
        @test ArbExtras._check_root_interval(sin, Arf(-1e-20), Arf(1e-20)) == (true, true)
        @test ArbExtras._check_root_interval(cos, Arf(1), Arf(2)) == (true, true)
        @test ArbExtras._check_root_interval(x -> sin(π * x / 2), Arf(1.5), Arf(2.5)) ==
              (true, true)
        @test ArbExtras._check_root_interval(
            x -> sin(π * x / 2),
            Arf(2 - 1e-10),
            Arf(2 + 1e-10),
        ) == (true, true)
        @test ArbExtras._check_root_interval(
            x -> (x - 1) * (x - 2) * (x - 3),
            Arf(1.9),
            Arf(2.1),
        ) == (true, true)

        # No zero in the interval
        @test ArbExtras._check_root_interval(sin, Arf(1e-10), Arf(2e-10)) == (false, false)
        @test ArbExtras._check_root_interval(sin, Arf(1), Arf(3)) == (false, false)
        @test ArbExtras._check_root_interval(cos, Arf(-1), Arf(1)) == (false, false)
        @test ArbExtras._check_root_interval(x -> sin(π * x / 2), Arf(0.5), Arf(1.5)) ==
              (false, false)
        @test ArbExtras._check_root_interval(
            x -> (x - 1) * (x - 2) * (x - 3),
            Arf(1.4),
            Arf(1.6),
        ) == (false, false)

        # Zero at endpoint
        @test ArbExtras._check_root_interval(sin, Arf(0), Arf(1)) == (true, false)
        @test ArbExtras._check_root_interval(sin, Arf(-1), Arf(0)) == (true, false)
        @test ArbExtras._check_root_interval(cospi, Arf(0), Arf(1 // 2)) == (true, false)

        # Multiple simple zeros in the interval
        @test ArbExtras._check_root_interval(sin, Arf(-1), Arf(4)) == (true, false)
        @test ArbExtras._check_root_interval(sin, Arf(-4), Arf(1)) == (true, false)
        @test ArbExtras._check_root_interval(cos, Arf(-2), Arf(2)) == (true, false)
        @test ArbExtras._check_root_interval(x -> sin(π * x / 2), Arf(1), Arf(2)) ==
              (true, false)
        @test ArbExtras._check_root_interval(x -> sin(π * x / 2), Arf(2), Arf(3)) ==
              (true, false)
        @test ArbExtras._check_root_interval(
            x -> (x - 1) * (x - 1 + eps()),
            Arf(1 - eps()),
            Arf(1 + eps()),
        ) == (true, false)

        # Zero of higher order
        @test ArbExtras._check_root_interval(x -> x^2, Arf(-eps()), Arf(eps())) ==
              (true, false)

        # check_unique = false
        @test ArbExtras._check_root_interval(sin, Arf(-1), Arf(1), check_unique = false) ==
              (true, false)
        @test ArbExtras._check_root_interval(cos, Arf(-1), Arf(1), check_unique = false) ==
              (false, false)
        @test ArbExtras._check_root_interval(sin, Arf(-1), Arf(4), check_unique = false) ==
              (true, false)
        @test ArbExtras._check_root_interval(
            x -> x^2,
            Arf(-eps()),
            Arf(eps()),
            check_unique = false,
        ) == (true, false)
    end

    @testset "ArbPoly" begin
        x = ArbPoly([1, 0])
        # Unique zero in the interval
        @test ArbExtras._check_root_interval(
            Arblib.fromroots(ArbPoly, [0]),
            Arf(-1),
            Arf(1),
        ) == (true, true)
        @test ArbExtras._check_root_interval(
            Arblib.fromroots(ArbPoly, [-2, 0, 2]),
            Arf(-1),
            Arf(1),
        ) == (true, true)
        @test ArbExtras._check_root_interval(
            Arblib.fromroots(ArbPoly, [0], [2im]),
            Arf(-1),
            Arf(1),
        ) == (true, true)

        # No zero in the interval
        @test ArbExtras._check_root_interval(ArbPoly(1), Arf(-1), Arf(1)) == (false, false)
        @test ArbExtras._check_root_interval(
            Arblib.fromroots(ArbPoly, [-2, 2]),
            Arf(-1),
            Arf(1),
        ) == (false, false)

        # Zero at endpoint
        @test ArbExtras._check_root_interval(ArbPoly([0, 1]), Arf(0), Arf(1)) ==
              (true, false)
        @test ArbExtras._check_root_interval(ArbPoly([0, 1]), Arf(-1), Arf(0)) ==
              (true, false)
        @test ArbExtras._check_root_interval(ArbPoly([-1, 0, 1]), Arf(0), Arf(1)) ==
              (true, false)

        # Multiple simple zeros in the interval
        @test ArbExtras._check_root_interval(
            Arblib.fromroots(ArbPoly, [-0.5, 0.5]),
            Arf(-1),
            Arf(1),
        ) == (true, false)
        @test ArbExtras._check_root_interval(
            Arblib.fromroots(ArbPoly, [-1, 1]),
            Arf(-1),
            Arf(1),
        ) == (true, false)

        # Zero of higher order
        @test ArbExtras._check_root_interval(ArbPoly(0), Arf(-1), Arf(1)) == (true, false)
        @test ArbExtras._check_root_interval(
            Arblib.fromroots(ArbPoly, [0, 0]),
            Arf(-1),
            Arf(1),
        ) == (true, false)

        # check_unique = false
        @test ArbExtras._check_root_interval(
            Arblib.fromroots(ArbPoly, [0]),
            Arf(-1),
            Arf(1),
            check_unique = false,
        ) == (true, false)
        @test ArbExtras._check_root_interval(
            ArbPoly(1),
            Arf(-1),
            Arf(1),
            check_unique = false,
        ) == (false, false)
        @test ArbExtras._check_root_interval(
            Arblib.fromroots(ArbPoly, [-0.5, 0.5]),
            Arf(-1),
            Arf(1),
            check_unique = false,
        ) == (true, false)
        @test ArbExtras._check_root_interval(
            ArbPoly(0),
            Arf(-1),
            Arf(1),
            check_unique = false,
        ) == (true, false)
    end
end

@testset "isolate_roots" begin
    @testset "sin(πx/2)" begin
        f = x -> sin(π * x / 2)
        a = -9.5
        b = 19.5
        roots = -8:2:18
        found, flags = ArbExtras.isolate_roots(f, Arf(a), Arf(b))

        # Check that all roots are accounted for
        found_balls = Arb.(found)
        @test all(any(contains(root), found_balls) for root in roots)

        # Check that all roots are reported as unique and that they
        # indeed are
        @test all(flags)
        @test all(isone(count(contains(root), found_balls)) for root in roots)
    end

    @testset "sin(1/x)" begin
        f = x -> sin(1 / x)
        a = 0
        b = 1
        found1, flags1 = ArbExtras.isolate_roots(f, Arf(a), Arf(b))
        found2, flags2 = ArbExtras.isolate_roots(f, Arf(a), Arf(b), depth = 15)

        # Check that the interval containing zero is included and that
        # it's smaller if we bisect more
        @test !flags1[1] && iszero(found1[1][1])
        @test !flags2[1] && iszero(found2[1][1])
        @test found1[1][2] > found2[1][2]

        # Check that we find isolate more zeros with higher depth
        @test sum(flags1) < sum(flags2)
    end

    @testset "thin interval" begin
        @test ArbExtras.isolate_roots(sin, Arf(0), Arf(0)) == ([(Arf(0), Arf(0))], [true])
        @test ArbExtras.isolate_roots(x -> sin(x) - sin(x), Arf(1), Arf(1)) ==
              ([(Arf(1), Arf(1))], [false])
        @test all(isempty.(ArbExtras.isolate_roots(cos, Arf(0), Arf(0))))
    end

    # Check that the derivative is not computed when check_unique is
    # false and that it splits the expected number of times
    @testset "check_unique = false" begin
        depth = 10
        found, flags = ArbExtras.isolate_roots(
            x -> Arb(0),
            Arf(0),
            Arf(1),
            check_unique = false,
            depth = depth,
        )
        @test length(found) == 2^(depth = 10 - 1)
        @test !any(flags)
        @test all(found[i][2] == found[i+1][1] for i = 1:length(found)-1)
    end

    @testset "ArbPoly" begin
        roots = Arb[-5, 0, ℯ, π, 10]
        p = Arblib.fromroots(ArbPoly, roots)
        a = -10
        b = 20
        found, flags = ArbExtras.isolate_roots(p, Arf(a), Arf(b), depth = 15)

        # Check that all roots are accounted for
        found_balls = Arb.(found)
        @test all(any(contains(root), found_balls) for root in roots)

        # Check that all roots are reported as unique and that they
        # indeed are
        @test all(flags)
        @test all(isone(count(contains(root), found_balls)) for root in roots)
    end

    @testset "ArbPoly - check_unique = false" begin
        depth = 10
        found, flags = ArbExtras.isolate_roots(
            ArbPoly(0),
            Arf(0),
            Arf(1),
            check_unique = false,
            depth = depth,
        )
        @test length(found) == 2^(depth = 10 - 1)
        @test !any(flags)
        @test all(found[i][2] == found[i+1][1] for i = 1:length(found)-1)
    end
end
