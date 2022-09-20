@testset "refine_root" begin
    # Define a number of problems for which we expect the method to
    # work without issues. The problems are defined as (f, enclosure,
    # answer, atol, rtol).
    problems = [
        # Problems with the root away from zero so that the relative
        # tolerance is relevant
        [
            (sin, cos, add_error(root, Mag(1)), root, 0, sqrt(eps(Arb))) for
            root in Arb(π) .* [-10:-1; 1:10]
        ]...,
        (
            x -> exp(x - 1) - 1,
            x -> exp(x - 1),
            setball(Arb, 1.5, 1),
            Arb(1),
            0,
            sqrt(eps(Arb)),
        ),
        (
            x -> exp(x + 1) - 1,
            x -> exp(x + 1),
            setball(Arb, -1.5, 1),
            Arb(-1),
            0,
            sqrt(eps(Arb)),
        ),

        # Problems where the root is at zero so that the absolute
        # tolerance is relevant
        (sin, cos, setball(Arb, 0.1, 1), Arb(0), sqrt(eps(Arb)), 0),
        (x -> exp(x) - 1, exp, setball(Arb, 0.1, 1), Arb(0), sqrt(eps(Arb)), 0),

        # ArbPoly
        (
            Arblib.fromroots(ArbPoly, [1, 2, 3]),
            Arblib.derivative(Arblib.fromroots(ArbPoly, [1, 2, 3])),
            setball(Arb, 2.001, 0.01),
            Arb(2),
            0,
            sqrt(eps(Arb)),
        ),
        (
            Arblib.fromroots(ArbPoly, [-Arb(π), Arb(π)]),
            Arblib.derivative(Arblib.fromroots(ArbPoly, [-Arb(π), Arb(π)])),
            add_error(Arb(π), Mag(0.1)),
            Arb(π),
            0,
            sqrt(eps(Arb)),
        ),
    ]

    for (f, df, enclosure, answer, atol, rtol) in problems
        root = ArbExtras.refine_root(f, enclosure; atol, rtol)
        aerror = 2radius(Arb, root)
        rerror = aerror / root
        @test Arblib.overlaps(root, answer)
        @test aerror < atol || rerror < rtol

        # Giving derivative explicitly
        root = ArbExtras.refine_root(f, enclosure; df, atol, rtol)
        aerror = 2radius(Arb, root)
        rerror = aerror / root
        @test Arblib.overlaps(root, answer)
        @test aerror < atol || rerror < rtol
    end

    # Test that it returns NaN in case it can't prove we have a root unless strict is true
    # Double root
    @test isnan(ArbExtras.refine_root(x -> cos(x) - 1, setball(Arb, 0, 1)))
    @test isequal(
        ArbExtras.refine_root(x -> cos(x) - 1, setball(Arb, 0, 1), strict = false),
        setball(Arb, 0, 1),
    )

    # Test that it stops correctly if the enclosure doesn't improve
    # and that it runs the minimum number of iterations still. We do
    # this by implementing a version of sin that counts the number of
    # times it's called.
    count = 0
    g = x -> begin
        count += 1
        return sin(x)
    end
    ArbExtras.refine_root(g, Arb(π), rtol = 0, min_iterations = 5)
    @test count == 10

    # Root on endpoint of interval
    # Enclosure which has 1 as an endpoint
    enclosure = setball(Arb, 1.125, Arblib.set_ui_2exp!(Mag(), unsigned(1), -3))
    @test isnan(ArbExtras.refine_root(sinpi, enclosure))
    root = ArbExtras.refine_root(sinpi, enclosure, strict = false)
    @test Arblib.overlaps(root, enclosure)
    @test contains(root, 1)

    # Test that it doesn't find a root if the root is just outside the
    # enclosure
    @test !contains(enclosure + Arb(1e-40), 1)
    @test isnan(ArbExtras.refine_root(sinpi, enclosure + Arb(1e-40)))

    # Check max_iterations

    # The function x -> x - 1 but with a weird implementation such
    # that for Arb it returns union(x, 1) - 1. Note that this still
    # gives a correct enclosure. This is done so that the convergence
    # for the Newton iterations is very slow.
    f(x::Arb) = begin
        return union(Arb(1), x) - 1
    end
    f(x) = x - 1

    @test radius(ArbExtras.refine_root(f, setball(Arb, 1.1, 1))) <
          radius(ArbExtras.refine_root(f, setball(Arb, 1.1, 1), atol = 1e-1))
    @test radius(ArbExtras.refine_root(f, setball(Arb, 1.1, 1), max_iterations = 40)) <
          radius(ArbExtras.refine_root(f, setball(Arb, 1.1, 1), max_iterations = 20))
end

@testset "refine_root_bisection" begin
    # Define a number of problems for which we expect the method to
    # work without issues. The problems are defined as (f, enclosure,
    # answer, atol, rtol).
    problems = [
        # Problems with the root away from zero so that the relative
        # tolerance is relevant
        [
            (sin, add_error(root, Mag(1)), root, 0, sqrt(sqrt(eps(Arb)))) for
            root in Arb(π) .* [-10:-1; 1:10]
        ]...,
        (x -> exp(x - 1) - 1, setball(Arb, 1.5, 1), Arb(1), 0, sqrt(sqrt(eps(Arb)))),
        (x -> exp(x + 1) - 1, setball(Arb, -1.5, 1), Arb(-1), 0, sqrt(sqrt(eps(Arb)))),

        # Problems where the root is at zero so that the absolute
        # tolerance is relevant
        (sin, setball(Arb, 0.1, 1), Arb(0), sqrt(sqrt(eps(Arb))), 0),
        (x -> exp(x) - 1, setball(Arb, 0.1, 1), Arb(0), sqrt(sqrt(eps(Arb))), 0),
        # Zero is exactly at the midpoint
        (sin, setball(Arb, 0, 1), Arb(0), sqrt(sqrt(eps(Arb))), 0),

        # ArbPoly
        (
            Arblib.fromroots(ArbPoly, [1, 2, 3]),
            setball(Arb, 2.001, 0.01),
            Arb(2),
            0,
            sqrt(sqrt(eps(Arb))),
        ),
        (
            Arblib.fromroots(ArbPoly, [-Arb(π), Arb(π)]),
            Arb((3, 4)),
            Arb(π),
            0,
            sqrt(sqrt(eps(Arb))),
        ),
    ]

    # Check that it works in the good cases
    for (f, enclosure, answer, atol, rtol) in problems
        root =
            Arb(ArbExtras.refine_root_bisection(f, getinterval(enclosure)...; atol, rtol))
        @test Arblib.overlaps(root, answer)
        @test ArbExtras.check_tolerance(root; atol, rtol)
    end

    # Test that it returns NaN in case the endpoints have the same sign unless strict = false
    @test all(isnan, ArbExtras.refine_root_bisection(x -> x^2, Arf(-1), Arf(1)))
    @test ArbExtras.refine_root_bisection(x -> x^2, Arf(-1), Arf(1), strict = false) ==
          (-1, 1)
    @test_logs (:warn, "sign of endpoints don't differ") match_mode = :any ArbExtras.refine_root_bisection(
        x -> x^2,
        Arf(-1),
        Arf(1),
        verbose = true,
    )

    # Test that it returns NaN in case the sign of the endpoints can't
    # be determined unless strict = false
    @test all(isnan, ArbExtras.refine_root_bisection(x -> x^2 - 1, Arf(-1), Arf(1)))
    @test ArbExtras.refine_root_bisection(x -> x^2 - 1, Arf(-1), Arf(1), strict = false) ==
          (-1, 1)
    @test_logs (:warn, "could not determine sign at right endpoint") match_mode = :any ArbExtras.refine_root_bisection(
        x -> x^2 - 1,
        Arf(-1),
        Arf(1),
        verbose = true,
    )

    # Check behaviour when sign of midpoint can't be determined
    @test_logs (:info, "could not determine sign at midpoint - shifting midpoint") match_mode =
        :any ArbExtras.refine_root_bisection(sin, Arf(-1), Arf(1), verbose = true)
    # Case when shifting midpoint also fails
    @test ArbExtras.refine_root_bisection(
        x -> ifelse(abs(x) < 1, zero(x), x),
        Arf(-1),
        Arf(1),
    ) == (-1, 1)
    @test_logs (:info, "could not determine sign at new midpoint") match_mode = :any ArbExtras.refine_root_bisection(
        x -> ifelse(abs(x) < 1, zero(x), x),
        Arf(-1),
        Arf(1),
        verbose = true,
    )

    # Check that it stops once the midpoint is equal to one of the
    # endpoints
    a, b = ArbExtras.refine_root_bisection(
        sin,
        Arf(3),
        Arf(4),
        rtol = 0,
        max_iterations = 1000,
    )
    @test (a + b) / 2 == a
    @test_logs (:info, "midpoint equal to endpoint, maximum precision reached") match_mode =
        :any ArbExtras.refine_root_bisection(
        sin,
        Arf(3),
        Arf(4),
        rtol = 0,
        max_iterations = 1000,
        verbose = true,
    )
end
