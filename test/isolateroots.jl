@testset "check_interval" begin
    # Once Arblib has better support for ArbSeries this should be set
    # to true
    check_unique = false
    @test_broken check_unique

    # Unique zero in the interval
    @test ArbExtras.check_interval(sin, Arf(-1), Arf(1); check_unique) == (true, false)
    @test ArbExtras.check_interval(sin, Arf(-1e-20), Arf(1e-20); check_unique) ==
          (true, false)
    @test ArbExtras.check_interval(cos, Arf(1), Arf(2); check_unique) == (true, false)
    @test ArbExtras.check_interval(x -> sin(π * x / 2), Arf(1.5), Arf(2.5); check_unique) ==
          (true, false)
    @test ArbExtras.check_interval(
        x -> sin(π * x / 2),
        Arf(2 - eps()),
        Arf(2 + eps());
        check_unique,
    ) == (true, false)
    @test ArbExtras.check_interval(
        x -> (x - 1) * (x - 2) * (x - 3),
        Arf(1.5),
        Arf(2.5);
        check_unique,
    ) == (true, false)

    # No zero in the interval
    @test ArbExtras.check_interval(sin, Arf(1e-10), Arf(2e-10); check_unique) ==
          (false, false)
    @test ArbExtras.check_interval(sin, Arf(1), Arf(3); check_unique) == (false, false)
    @test ArbExtras.check_interval(cos, Arf(-1), Arf(1); check_unique) == (false, false)
    @test ArbExtras.check_interval(x -> sin(π * x / 2), Arf(0.5), Arf(1.5); check_unique) ==
          (false, false)
    @test ArbExtras.check_interval(
        x -> (x - 1) * (x - 2) * (x - 3),
        Arf(1.4),
        Arf(1.6);
        check_unique,
    ) == (false, false)

    # Multiple simple zeros in the interval
    @test ArbExtras.check_interval(sin, Arf(-1), Arf(4); check_unique) == (true, false)
    @test ArbExtras.check_interval(sin, Arf(-4), Arf(1); check_unique) == (true, false)
    @test ArbExtras.check_interval(cos, Arf(-2), Arf(2); check_unique) == (true, false)
    @test ArbExtras.check_interval(x -> sin(π * x / 2), Arf(1), Arf(2); check_unique) ==
          (true, false)
    @test ArbExtras.check_interval(x -> sin(π * x / 2), Arf(2), Arf(3); check_unique) ==
          (true, false)
    @test ArbExtras.check_interval(
        x -> (x - 1) * (x - 1 + eps()),
        Arf(1 - eps()),
        Arf(1 + eps()),
        ;
        check_unique,
    ) == (true, false)
    @test ArbExtras.check_interval(x -> x^2, Arf(-eps()), Arf(eps())) == (true, false)
end

@testset "isolateroots" begin
    check_unique = false
    @test_broken check_unique

    f = x -> sin(π * x / 2)
    a = -10
    b = 20
    found, flags = ArbExtras.isolateroots(f, Arf(-10), Arf(20); check_unique)

    # Check that all roots are accounted for
    found_balls = Arb.(found)
    @test all(any(contains.(found_balls, i)) for i = a:2:b)

    # Check that all reported unique roots are unique by check that
    # their radius is smaller than 1 (the distance between the roots
    # being 2)
    # TODO: The result could be okay even though this fails
    for i in eachindex(found)
        if flags[i]
            @test Arf.(Arblib.radref(found_balls[i])) < 1
        end
    end

    # It should be able to successfully find all roots in this case
    @test_broken sum(flags) == length(a:2:b)
end
