demo_problems_definitions = [
    (x -> sin(x) + sin(10 // 3 * x), (2.7, 7.5), 5.145735, 6.2173, true),
    (x -> -sum(k * sin((k + 1) * x + k) for k = 1:6), (-10, 10), -6.84128, -1.0991, true),
    (x -> -(16x^2 - 24x + 5) * exp(-x), (1.9, 3.9), 2.868034, :left, false),
    (x -> -x * (1.4 - 3x) * sin(18x), (0, 1.2), 0.969171, 1.1416, true),
    (x -> -(x + sin(x)) * exp(-x^2), (-10, 10), 0.67956, -0.6790, true),
    (x -> sin(x) + sin(10 // 3 * x) + log(x) - 0.84x + 3, (2.7, 7.5), 5.19978, :left, true),
    (x -> -sum(k * cos((k + 1) * x + k) for k = 1:6), (-10, 10), -7.1095, 4.9134, true),
    (x -> sin(x) + sin(2 // 3 * x), (3.1, 20.4), 17.039, :right, true),
    (x -> -x * sin(x), (0, 10), 7.9787, :right, true),
    (x -> 2cos(x) + cos(2x), (-π / 2, 2π), 2.09439, 2π, true), # Handle
    (x -> sin(x)^3 + cos(x)^3, (0, 2π), π, :left, true),
    (x -> -x^(2 // 3) - (1 - x^2)^(1 // 3), (0.001, 0.99), 1 / sqrt(2), :left, false),
    (x -> -exp(-x) * sinpi(2x), (0, 4), 0.224885, 0.7248, true),
    (x -> (x^2 - 5x + 6) / (x^2 + 1), (-5, 5), 2.41422, -0.4145, true),
    (x -> -(x - sin(x)) * exp(-x^2), (-10, 10), 1.195137, -1.1951, true),
    (x -> x * sin(x) + x * cos(2x), (0, 10), 4.79507, 9.2039, true),
    (x -> exp(-3x) - sin(x)^3, (0, 20), 9π / 2, 4.7124, true),
]

"""
    demo_problem(i, thin = false)

Return a demo problem for testing `extrema_enclosure` and
`extrema_series`. Returns `(f, (a, b), fmin, fmax, fabsmin, fabsmax)`
where `f` is the function, `(a, b)` gives the interval to use, `fmin`
and `fmax` are the minimum and maximum value on the interval and
`fabsmin` and `fabsmax` are the minimum and maximum value of the
absolute value of `f`.

If `thin = true` then don't return the full interval `[a, b]` but a
thinner one around the location of the minimum.

The problems are taken from
http://infinity77.net/global_optimization/test_functions_1d.html. The
minima given there are only approximate and not even correct in all
cases. To get high precision values we take an a priori approximation
for the location of the minimum and then we rigorously find a zero of
the derivative close to that approximation.
"""
function demo_problem(i::Integer, thin = false)
    f, (a, b), xmin, xmax, crosses = demo_problems_definitions[i]

    a, b = convert(Arf, a), convert(Arf, b)

    if xmin == :left
        xmin = Arb(a)
    elseif xmin == :right
        xmin = Arb(b)
    else
        xmin = add_error(convert(Arb, xmin), Mag(1e-3))
    end

    if xmax == :left
        xmax = Arb(a)
    elseif xmax == :right
        xmax = Arb(b)
    else
        xmax = add_error(convert(Arb, xmax), Mag(1e-3))
    end


    df(x::Arb) = f(ArbSeries((x, 1), degree = 1))[1]
    df(x::ArbSeries) = Arblib.derivative(f(ArbSeries(x, degree = Arblib.degree(x) + 1)))

    if !iszero(Arblib.radref(xmin))
        roots, flags = ArbExtras.isolate_roots(df, Arblib.getinterval(Arf, xmin)...)
        @assert only(flags)
        root_min = ArbExtras.refine_root(df, Arb(only(roots)))
        xmin = Arblib.intersect(Arb((a, b)), root_min)
    end

    if !iszero(Arblib.radref(xmax))
        roots, flags = ArbExtras.isolate_roots(df, Arblib.getinterval(Arf, xmax)...)
        @assert only(flags)
        root_max = ArbExtras.refine_root(df, Arb(only(roots)))
        xmax = Arblib.intersect(Arb((a, b)), root_max)
    end

    if thin
        # Thin interval around roots
        interval_min = add_error(xmin, Mag(1e-2))
        interval_max = add_error(xmax, Mag(1e-2))
        # Intersect with interval
        c_min, d_min = getinterval(Arf, intersect(Arb((a, b)), interval_min))
        c_max, d_max = getinterval(Arf, intersect(Arb((a, b)), interval_max))

        return f, (c_min, d_min), f(xmin)
    else
        absmin = crosses ? zero(Arb) : min(abs(f(xmin)), abs(f(xmax)))
        absmax = max(abs(f(xmin)), abs(f(xmax)))
        return f, (a, b), f(xmin), f(xmax), absmin, absmax
    end
end
