"""
    integrate_gauss_legendre(f, a::Arb, b::Arb)

Compute an enclosure of the integral of `f` on the interval `[a, b]`
using a Gauss-Legendre quadrature of order 2.

The approximation is given by
```
(b - a) / 2 * (f(x₁) + f(x₂))
```
where
```
x₁ = (b - a) / 2 * sqrt(3) / 3 + (b + a) / 2
x₂ = -(b - a) / 2 * sqrt(3) / 3 + (b + a) / 2
```
are quadrature nodes. The remainder term is enclosed by
```
(b - a)^5 / 4320 * d4f
```
where `d4f` is an enclosure of the fourth derivative of `f` on the
interval `[a, b]`.

If it happens to be that `d4f` is not finite, for example because `f`
is not differentiable on the interval, it instead computes a naive
enclosure of the integral using the enclosure of `f` evaluated on the
interval `[a, b]`.

The function `f` should support evaluation on both `Arb` and
`ArbSeries` and should return an enclosure of the result in both
cases. For `ArbSeries` it is enough that it supports evaluation up to
degree 4 and only the zeroth and fourth order term are ever used.

- **TODO:** Optimize for performance
- **TODO:** Add option for computing integral of absolute value?
- **TODO:** Allow giving an separate function for computing derivative?
"""
function integrate_gauss_legendre(f, a::Arb, b::Arb)
    bma = b - a

    # x_series = ArbSeries((Arblib.union(a, b), 1), degree = 4)
    x_series = ArbSeries(degree = 4, prec = precision(bma))
    # We need to do it in this order so that the degree of the
    # polynomial is correct when we set the constant coefficient.
    x_series[1] = 1
    Arblib.union!(Arblib.ref(x_series, 0), a, b)

    f_series = f(x_series)

    # If the fourth derivative of f is not finite compute an enclosure
    # of the integral using the enclosure of f on the interval [a, b].
    isfinite(Arblib.ref(f_series, 4)) ||
        return Arblib.mul!(bma, bma, Arblib.ref(f_series, 0))

    # remainder = (b - a)^5 / 4320 * f_series[4] * factorial(4)
    remainder = Arblib.pow!(zero(bma), bma, UInt(5))
    Arblib.mul!(remainder, remainder, Arblib.ref(f_series, 4))
    Arblib.div!(remainder, remainder, 180)

    # v = (b - a) / 2 * sqrt(Arb(3)) / 3
    v = Arblib.sqrt!(zero(bma), UInt(3))
    Arblib.mul!(v, v, bma)
    Arblib.div!(v, v, 6)

    # Set x₁ = (b + a) / 2 to begin with
    x₁ = b + a
    Arblib.mul_2exp!(x₁, x₁, -1)

    # x₂ = (b + a) / 2 - v = x₁ - v
    x₂ = Arblib.sub!(zero(bma), x₁, v)

    # x₁ = (b + a) / 2 - v = x₁ - v
    Arblib.add!(x₁, x₁, v)

    # approximation = (b - a) / 2 * (f(x₁) + f(x₂))
    approximation = f(x₁)
    Arblib.add!(approximation, approximation, f(x₂))
    Arblib.mul!(approximation, approximation, bma)
    Arblib.mul_2exp!(approximation, approximation, -1)

    # Return approximation + remainder
    return Arblib.add!(approximation, approximation, remainder)
end

"""
    integrate(f, a::Arb, b::Arb; atol, rtol, depth_start, maxevals, depth, verbose)

Compute an enclosure of the integral of `f on the interval `[a, b]`.

It uses a Guass-Legendre quadrature of order 2 through
[`integrate_gauss_legendre`](@ref) together with bisection of the
interval. On each subinterval the integral is computed using
[`integrate_gauss_legendre`](@ref) and it checks if the result
satisfies the required tolerance, if it doesn't the interval is
bisected. For more details about the integration method see
[`integrate_gauss_legendre`](@ref).

Notice that the tolerance is checked on each **subinterval** and not
on the interval as a whole. This means that even if it finishes before
reaching the maximum depth or the maximum number of evaluations the
result as a whole might not satisfy the given tolerance.

The argument `depth_start` bisect the interval using
[`bisect_interval_recursive`](@ref) before starting to compute the
integral. This can be useful if it is known beforehand that a certain
number of bisections will be necessary before the enclosures get good
enough. It defaults to `0` which corresponds to not bisecting the
interval at all before starting.

The arguments `maxevals` and `depth` can be used to limit the number
of function evaluations and the number of bisections of the interval
respectively.

If `verbose = true` then output information about the process.
"""
function integrate(
    f,
    a::Arb,
    b::Arb;
    atol = 0,
    rtol = sqrt(eps(one(a))),
    depth_start::Integer = 0,
    maxevals::Integer = 1000,
    depth::Integer = 20,
    verbose = false,
)
    isfinite(a) && isfinite(b) || return Arblib.indeterminate!(zero(a))

    intervals = bisect_interval_recursive(a, b, depth_start)

    integral = zero(a)

    verbose && @info "iteration: $(lpad(0, 2)), " *
          " starting intervals: $(lpad(length(intervals), 3))"

    iterations = 0
    evals = 0
    while !isempty(intervals)
        iterations += 1
        evals += length(intervals)

        # If we have reached the maximum number of iterations or
        # evaluations then we don't split any further.
        finished = iterations >= depth || evals >= maxevals

        # The number of intervals that had to be split further. The
        # only reason we store this in a separate variable is to
        # correctly count the number during the last iteration when
        # finished = true.
        remaining_intervals = 0

        to_split = falses(length(intervals))
        for (i, (a, b)) in enumerate(intervals)
            integral_part = integrate_gauss_legendre(f, a, b)

            if check_tolerance(integral_part; atol, rtol)
                integral += integral_part
            elseif finished
                remaining_intervals += 1
                integral += integral_part
            else
                remaining_intervals += 1
                to_split[i] = true
            end
        end
        intervals = bisect_intervals(intervals, to_split)

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(remaining_intervals, 3)), " *
              "current integral: $integral"
    end

    return integral
end
