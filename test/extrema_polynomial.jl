@testset "_extrema_polynomial_low_degree" begin
    # Degree 0
    for p0 in (Arb(0), Arb(-1), Arb(1), Arb((-1, 1)))
        @test isequal(
            ArbExtras._extrema_polynomial_low_degree(ArbPoly(p0), Arf(-1), Arf(1)),
            (p0, p0),
        )

        @test isequal(
            ArbExtras._minimum_polynomial_low_degree(ArbPoly(p0), Arf(-1), Arf(1)),
            p0,
        )

        @test isequal(
            ArbExtras._maximum_polynomial_low_degree(ArbPoly(p0), Arf(-1), Arf(1)),
            p0,
        )
    end

    # Degree 1
    let p = ArbPoly((1, 3))
        # abs_value = false
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(1)) == (-2, 4)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(0), Arf(1)) == (1, 4)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-2), Arf(-1)) == (-5, -2)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(1)) == (-4, 2)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(0), Arf(1)) == (-4, -1)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-2), Arf(-1)) == (2, 5)

        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(1)) == -2
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(0), Arf(1)) == 1
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-2), Arf(-1)) == -5
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(1)) == -4
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(0), Arf(1)) == -4
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-2), Arf(-1)) == 2

        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(1)) == 4
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(0), Arf(1)) == 4
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-2), Arf(-1)) == -2
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(1)) == 2
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(0), Arf(1)) == -1
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-2), Arf(-1)) == 5

        # abs_value = true
        abs_value = true
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value) ==
              (0, 4)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(0), Arf(1); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(0), Arf(1); abs_value) ==
              (1, 4)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-2), Arf(-1); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-2), Arf(-1); abs_value) ==
              (2, 5)

        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(0), Arf(1); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(0), Arf(1); abs_value) ==
              1
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-2), Arf(-1); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-2), Arf(-1); abs_value) ==
              2

        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value) ==
              4
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(0), Arf(1); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(0), Arf(1); abs_value) ==
              4
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-2), Arf(-1); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-2), Arf(-1); abs_value) ==
              5
    end

    let p = ArbPoly((1, Arb((-1, 2))))
        # abs_value = false
        @test all(
            contains.(
                ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(1)),
                (Arb((-1, 2)), Arb((0, 3))),
            ),
        )
        @test contains(
            ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(1)),
            Arb((-1, 2)),
        )
        @test contains(
            ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(1)),
            Arb((0, 3)),
        )
        @test all(
            contains.(
                ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(1)),
                (Arb((-3, 0)), Arb((-2, 1))),
            ),
        )
        @test contains(
            ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(1)),
            Arb((-3, 0)),
        )
        @test contains(
            ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(1)),
            Arb((-2, 1)),
        )

        # abs_value = true
        abs_value = true
        @test all(
            contains.(
                ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value),
                (Arb((0, 2)), Arb((0, 3))),
            ),
        )
        @test contains(
            ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value),
            Arb((0, 2)),
        )
        @test contains(
            ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value),
            Arb((0, 3)),
        )
        @test all(
            contains.(
                ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value),
                (Arb((0, 2)), Arb((0, 3))),
            ),
        )
        @test contains(
            ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value),
            Arb((0, 2)),
        )
        @test contains(
            ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value),
            Arb((0, 3)),
        )
    end

    # Degree 2
    let p = ArbPoly((-3, -1, 2))
        # abs_value = false
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(-2)) == (7, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(-1)) == (0, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(0)) == (-3, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(1)) == (-25 // 8, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(2)) == (-25 // 8, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(3)) == (-25 // 8, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(0)) == (-3, 0)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(1)) == (-25 // 8, 0)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(2)) == (-25 // 8, 3)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(3)) == (-25 // 8, 12)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(0), Arf(1)) == (-25 // 8, -2)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(0), Arf(2)) == (-25 // 8, 3)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(1), Arf(2)) == (-2, 3)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(1), Arf(3)) == (-2, 12)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(-2)) == (-18, -7)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(-1)) == (-18, -0)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(0)) == (-18, 3)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(1)) ==
              (-18, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(2)) ==
              (-18, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(3)) ==
              (-18, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(0)) == (-0, 3)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(1)) == (-0, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(2)) == (-3, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(3)) ==
              (-12, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(0), Arf(1)) == (2, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(0), Arf(2)) == (-3, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(1), Arf(2)) == (-3, 2)
        @test ArbExtras._extrema_polynomial_low_degree(-p, Arf(1), Arf(3)) == (-12, 2)

        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(-2)) == 7
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(-1)) == 0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(0)) == -3
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(1)) == -25 // 8
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(2)) == -25 // 8
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(3)) == -25 // 8
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(0)) == -3
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(1)) == -25 // 8
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(2)) == -25 // 8
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(3)) == -25 // 8
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(0), Arf(1)) == -25 // 8
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(0), Arf(2)) == -25 // 8
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(1), Arf(2)) == -2
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(1), Arf(3)) == -2
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(-2)) == -18
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(-1)) == -18
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(0)) == -18
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(1)) == -18
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(2)) == -18
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(3)) == -18
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(0)) == -0
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(1)) == -0
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(2)) == -3
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(3)) == -12
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(0), Arf(1)) == 2
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(0), Arf(2)) == -3
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(1), Arf(2)) == -3
        @test ArbExtras._minimum_polynomial_low_degree(-p, Arf(1), Arf(3)) == -12

        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(-2)) == 18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(-1)) == 18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(0)) == 18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(1)) == 18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(2)) == 18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(3)) == 18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(0)) == 0
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(1)) == 0
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(2)) == 3
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(3)) == 12
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(0), Arf(1)) == -2
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(0), Arf(2)) == 3
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(1), Arf(2)) == 3
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(1), Arf(3)) == 12
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(-2)) == -7
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(-1)) == -0
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(0)) == 3
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(1)) == 25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(2)) == 25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(3)) == 25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(0)) == 3
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(1)) == 25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(2)) == 25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(3)) == 25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(0), Arf(1)) == 25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(0), Arf(2)) == 25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(1), Arf(2)) == 2
        @test ArbExtras._maximum_polynomial_low_degree(-p, Arf(1), Arf(3)) == 2

        # abs_value = true
        abs_value = true
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(-2); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(-2); abs_value) ==
              (7, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(-1); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(-1); abs_value) ==
              (0, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(0); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(0); abs_value) ==
              (0, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(1); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(1); abs_value) ==
              (0, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(2); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(2); abs_value) ==
              (0, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-3), Arf(3); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-3), Arf(3); abs_value) ==
              (0, 18)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(0); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(0); abs_value) ==
              (0, 3)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value) ==
              (0, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(2); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(2); abs_value) ==
              (0, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(3); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(3); abs_value) ==
              (0, 12)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(0), Arf(1); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(0), Arf(1); abs_value) ==
              (2, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(0), Arf(2); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(0), Arf(2); abs_value) ==
              (0, 25 // 8)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(1), Arf(2); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(1), Arf(2); abs_value) ==
              (0, 3)
        @test ArbExtras._extrema_polynomial_low_degree(p, Arf(1), Arf(3); abs_value) ==
              ArbExtras._extrema_polynomial_low_degree(-p, Arf(1), Arf(3); abs_value) ==
              (0, 12)

        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(-2); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(-2); abs_value) ==
              7
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(-1); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(-1); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(0); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(0); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(1); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(1); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(2); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(2); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-3), Arf(3); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-3), Arf(3); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(0); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(0); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(2); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(2); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(3); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(3); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(0), Arf(1); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(0), Arf(1); abs_value) ==
              2
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(0), Arf(2); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(0), Arf(2); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(1), Arf(2); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(1), Arf(2); abs_value) ==
              0
        @test ArbExtras._minimum_polynomial_low_degree(p, Arf(1), Arf(3); abs_value) ==
              ArbExtras._minimum_polynomial_low_degree(-p, Arf(1), Arf(3); abs_value) ==
              0

        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(-2); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(-2); abs_value) ==
              18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(-1); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(-1); abs_value) ==
              18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(0); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(0); abs_value) ==
              18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(1); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(1); abs_value) ==
              18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(2); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(2); abs_value) ==
              18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-3), Arf(3); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-3), Arf(3); abs_value) ==
              18
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(0); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(0); abs_value) ==
              3
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value) ==
              25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(2); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(2); abs_value) ==
              25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(3); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(3); abs_value) ==
              12
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(0), Arf(1); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(0), Arf(1); abs_value) ==
              25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(0), Arf(2); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(0), Arf(2); abs_value) ==
              25 // 8
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(1), Arf(2); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(1), Arf(2); abs_value) ==
              3
        @test ArbExtras._maximum_polynomial_low_degree(p, Arf(1), Arf(3); abs_value) ==
              ArbExtras._maximum_polynomial_low_degree(-p, Arf(1), Arf(3); abs_value) ==
              12
    end

    let p = ArbPoly((-3, -1, Arb((-1, 2))))
        # abs_value = false
        @test all(
            contains.(
                ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(1)),
                (Arb((-5, -2)), Arb((-3, 0))),
            ),
        )
        @test contains(
            ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(1)),
            Arb((-5, -2)),
        )
        @test contains(
            ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(1)),
            Arb((-3, 0)),
        )
        @test all(
            contains.(
                ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(1)),
                (Arb((0, 3)), Arb((2, 5))),
            ),
        )
        @test contains(
            ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(1)),
            Arb((0, 3)),
        )
        @test contains(
            ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(1)),
            Arb((2, 5)),
        )

        # abs_value = true
        abs_value = true
        @test all(
            contains.(
                ArbExtras._extrema_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value),
                (Arb((0, 3)), Arb((2, 5))),
            ),
        )
        @test contains(
            ArbExtras._minimum_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value),
            Arb((0, 3)),
        )
        @test contains(
            ArbExtras._maximum_polynomial_low_degree(p, Arf(-1), Arf(1); abs_value),
            Arb((2, 5)),
        )
        @test all(
            contains.(
                ArbExtras._extrema_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value),
                (Arb((0, 3)), Arb((2, 5))),
            ),
        )
        @test contains(
            ArbExtras._minimum_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value),
            Arb((0, 3)),
        )
        @test contains(
            ArbExtras._maximum_polynomial_low_degree(-p, Arf(-1), Arf(1); abs_value),
            Arb((2, 5)),
        )
    end

    # Higher order degree
    @test_throws ArgumentError ArbExtras._extrema_polynomial_low_degree(
        ArbPoly((0, 0, 0, 1)),
        Arf(-1),
        Arf(1),
    )
end

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

    # Wide polynomials
    @test_logs (:info, "could not evaluate p to any meaningful precision") match_mode = :any ArbExtras.extrema_polynomial(
        ArbPoly(add_error.((Arb(1), Arb(1), Arb(1), Arb(1)), Mag(100))),
        Arf(-1),
        Arf(1),
        verbose = true,
    )
    @test_logs (:info, "could not evaluate p to any meaningful precision") match_mode = :any ArbExtras.minimum_polynomial(
        ArbPoly(add_error.((Arb(1), Arb(1), Arb(1), Arb(1)), Mag(100))),
        Arf(-1),
        Arf(1),
        verbose = true,
    )
    @test_logs (:info, "could not evaluate p to any meaningful precision") match_mode = :any ArbExtras.maximum_polynomial(
        ArbPoly(add_error.((Arb(1), Arb(1), Arb(1), Arb(1)), Mag(100))),
        Arf(-1),
        Arf(1),
        verbose = true,
    )

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
