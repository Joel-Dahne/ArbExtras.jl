@testset "refine_root" begin
    # Define a number of problems for which we expect the method to
    # work without issues. The problems are defined as (f, enclosure,
    # answer, atol, rtol).
    problems = [
        # Problems with the root away from zero so that the relative
        # tolerance is relevant
        [
            (sin, Arblib.add_error!(root, Mag(1)), root, 0, sqrt(eps(Arb))) for
            root in Arb(π) .* [-10:-1; 1:10]
        ]...,
        (
            x -> exp(x - 1) - 1,
            Arblib.add_error!(Arb(1.5), Mag(1)),
            Arb(1),
            0,
            sqrt(eps(Arb)),
        ),
        (
            x -> exp(x + 1) - 1,
            Arblib.add_error!(Arb(-1.5), Mag(1)),
            Arb(-1),
            0,
            sqrt(eps(Arb)),
        ),

        # Problems where the root is at zero so that the absolute
        # tolerance is relevant
        (sin, Arblib.add_error!(Arb(0.1), Mag(1)), Arb(0), sqrt(eps(Arb)), 0),
        (x -> exp(x) - 1, Arblib.add_error!(Arb(0.1), Mag(1)), Arb(0), sqrt(eps(Arb)), 0),

        # ArbPoly
        (
            Arblib.fromroots(ArbPoly, [1, 2, 3]),
            Arblib.add_error!(Arb(2.001), Mag(0.01)),
            Arb(2),
            0,
            sqrt(eps(Arb)),
        ),
        (
            Arblib.fromroots(ArbPoly, [-Arb(π), Arb(π)]),
            Arblib.add_error!(Arb(π), Mag(0.1)),
            Arb(π),
            0,
            sqrt(eps(Arb)),
        ),
    ]

    for (f, enclosure, answer, atol, rtol) in problems
        root = ArbExtras.refine_root(f, enclosure; atol, rtol)
        aerror = 2Arblib.radius(Arb, root)
        rerror = aerror / root
        @test Arblib.overlaps(root, answer)
        @test aerror < atol || rerror < rtol
    end

    # Test that it returns NaN in case it can't prove we have a root unless strict is true
    # Double root
    @test isnan(ArbExtras.refine_root(x -> cos(x) - 1, Arblib.add_error!(Arb(0), Mag(1))))
    @test isequal(
        ArbExtras.refine_root(
            x -> cos(x) - 1,
            Arblib.add_error!(Arb(0), Mag(1)),
            strict = false,
        ),
        Arblib.add_error!(Arb(0), Mag(1)),
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
    enclosure = Arblib.add_error!(Arb(1.125), Arblib.set_ui_2exp!(Mag(), unsigned(1), -3))
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

    @test Arblib.radius(ArbExtras.refine_root(f, Arblib.add_error!(Arb(1.1), Mag(1)))) <
          Arblib.radius(
        ArbExtras.refine_root(f, Arblib.add_error!(Arb(1.1), Mag(1)), atol = 1e-1),
    )
    @test Arblib.radius(
        ArbExtras.refine_root(f, Arblib.add_error!(Arb(1.1), Mag(1)), max_iterations = 40),
    ) < Arblib.radius(
        ArbExtras.refine_root(f, Arblib.add_error!(Arb(1.1), Mag(1)), max_iterations = 20),
    )
end
