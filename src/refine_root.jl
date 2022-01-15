"""
    refine_root(f, root::Arb; atol, rtol, min_iterations, max_iterations, strict)
    refine_root(f::ArbPoly, root::Arb; atol, rtol, min_iterations, max_iterations, strict)

Refine the given root of the function `f`.

Given a ball `root` containing a unique root of the function `f` it
refines the enclosure of the root using successive interval Newton
iterations. The method can only handle simple roots contained in the
interior of the enclosure and to work the derivative must be non-zero
on the whole enclosure.

The function `f` can be either an `ArbPoly` or a regular function. If
it's a regular function then the derivative is computed using
`ArbSeries`.

At each iteration it checks if the required tolerance is met according
to [`check_tolerance`](@ref). The default tolerances are `atol = 0`
and `rtol = 4eps(one(root))`.

If the absolute error does not improve by at least a factor 1.5
between two iterations then it stops early since it's unlikely that
further iterations would improve the result much. This can for example
happen when the function is computed with too low precision. This
check is only done if `min_iterations` have been performed, this is to
avoid stopping early due to slow convergence in the beginning. To
avoid this check entirely `min_iterations` can be put to the same as
`max_iterations`.

If `strict = true` then only return an enclosure which has been proved
to contain a unique root contained in the original enclosure, if it
could not prove that this is the case it returns `NaN`. If `strict =
false` it will return an enclosure in all cases, this enclosure is
only valid if `root` indeed contains a unique root and the returned
enclosure could just be the starting enclosure `root` itself. This is
useful if you know that `root` contains a unique root and don't want
to worry about this method failing.

If any iteration gives `NaN` as a result then it returns the enclosure
from the previous iteration, possibly just the starting enclosure.
This is unless `strict = true` and it has not been able to prove that
there is a root, in which case `NaN` is returned.

If `verbose = true` then print the enclosure at each iteration and
some more information in the end.

In rare cases when `strict = false` and the original enclosure does in
fact not contain a root the method could converge to a potential root
just outside the original enclosure. The reason this might happen is
that when the intersection of two balls is computed the enclosure is
typically slightly larger than the true intersection. This is not an
issue if `strict = false` nor if the original enclosure does contain a
root.
"""
function refine_root(
    f,
    root::Arb;
    atol = 0,
    rtol = 4eps(one(root)),
    min_iterations = 1,
    max_iterations = 20,
    strict = true,
    verbose = false,
)
    original_root = root
    root = copy(root)
    mid = Arblib.midpoint(Arb, root)
    if f isa ArbPoly
        df = Arblib.derivative(f)
    else
        series = ArbSeries((root, 1))
    end

    verbose && @info "enclosure: $root"

    error_previous = radius(root)
    isproved = false
    for i = 1:max_iterations
        # Compute new enclosure
        y = f(mid)
        if f isa ArbPoly
            dy = df(root)
        else
            dy = f(series)[1]
        end
        new_root = mid - y / dy

        # Note that since Arblib.intersection! only returns an
        # enclosure of the result it's not enough to check that the
        # new enclosure is contained in the interior of previous
        # enclosure. It could happen that the previous enclosure
        # contains numbers not contained in the original enclosure and
        # in that case we are not guaranteed that the root is
        # contained in the original enclosure, only in the previous
        # one.
        isproved |= Arblib.contains_interior(original_root, new_root)

        if isnan(new_root) || !Arblib.overlaps(root, new_root)
            break
        end

        Arblib.intersection!(root, root, new_root)
        verbose && @info "enclosure: $root"

        # If the result satisfies the required tolerance - break
        if check_tolerance(root; atol, rtol)
            verbose && @info "tolerance satisfied"
            break
        end

        # If the result did not improve compared to the last iteration
        # and we have performed the minimum number of iterations -
        # break
        error = radius(root)
        if i >= min_iterations && 1.5error > error_previous
            verbose &&
                @info "diameter only improved from $error_previous to $error - stopping early"
            break
        end
        error_previous = error

        # Update values
        Arblib.set!(mid, Arblib.midref(root))
        if !(f isa ArbPoly)
            series[0] = root
        end
    end

    if strict && !isproved
        verbose && @warn "could not prove root"
        return Arb(NaN, prec = precision(root))
    else
        return root
    end
end

"""
    refine_root_bisection(f, root::Arb; atol, rtol, max_iterations, strict, verbose)

Refine the root of the function `f` on the interval `[a, b]` using
bisection.

The function `f` is assumed to be continuous on the interval `[a, b]`
and `f(a)` and `f(b)` should have different signs.

It iteratively bisects the interval and determines which subinterval
to keep by checking the sign of `f` at the midpoint. It stops either
when the required tolerance is met according to
[`check_tolerance`](@ref), it has reached the maximum number of
iterations, it can't determine the sign of the midpoint or the
midpoint is equal to one of the endpoints (meaning that we have
bisected as much as possible at the given precision).

If the zero lies exactly at or very close to one of the bisection
point it can fail to determine the sign there. We avoid this by in
that case instead trying to bisect at a point between the left
endpoint of the interval and the midpoint, if that also fails it
stops.

When possible [`refine_root`](@ref) should be used instead since it
converges much faster. This method is useful when the derivative is
not available, i.e. evaluation with `ArbSeries` is not possible, or as
a preprocessor for [`refine_root`](@ref) if the original enclosure of
the root gives an enclosure on the derivative which contains zero.

If `strict = true` return `NaN` if the signs of `f(a)` and `f(b)`
don't differ. In this case a finite result proves that there is a root
of `f` on the interval, assuming that `f` is continuous there, but it
says nothing about the uniqueness.

If `verbose = true` then print the enclosure at each iteration and
some more information in the end.

- **IMPROVE:** Reduce the number of allocations.
- **IMPROVE:** Handle the case when the sign at the midpoint can't be
  determined better by using a perturbation of the midpoint.
"""
function refine_root_bisection(
    f,
    a::Arf,
    b::Arf;
    atol = 0,
    rtol = sqrt(eps(one(a))),
    max_iterations = precision(a) ÷ 2,
    strict = true,
    verbose = false,
)
    check_interval(a, b)

    sign_a, sign_b = Arblib.sgn_nonzero(f(Arb(a))), Arblib.sgn_nonzero(f(Arb(b)))

    if sign_a * sign_b > 0
        verbose && @warn "sign of endpoints don't differ"
        return strict ? (Arf(NaN), Arf(NaN)) : (a, b)
    end

    if sign_a == 0 || sign_b == 0
        verbose && sign_a == 0 && @warn "could not determine sign at left endpoint"
        verbose && sign_b == 0 && @warn "could not determine sign at right endpoint"
        return strict ? (Arf(NaN), Arf(NaN)) : (a, b)
    end

    root = Arb((a, b)) # Compute root to be able to check tolerance

    for iteration = 1:max_iterations
        if check_tolerance(root; atol, rtol)
            verbose && @info "tolerance satisfied"
            break
        end

        mid = Arblib.mul_2exp!(zero(Arf), a + b, -1)

        if a == mid || b == mid
            verbose && @info "midpoint equal to endpoint, maximum precision reached"
            break
        end

        sign_mid = Arblib.sgn_nonzero(f(Arb(mid)))

        if sign_mid == 0
            verbose && @info "could not determine sign at midpoint - shifting midpoint"

            # Try with a different midpoint
            mid = Arblib.mul_2exp!(zero(Arf), a + mid, -1)
            sign_mid = Arblib.sgn_nonzero(f(Arb(mid)))

            # If this also fails we stop
            if sign_mid == 0
                verbose && @info "could not determine sign at new midpoint"
                break
            end
        end

        # We always have sign_mid != 0 here

        if sign_mid == sign_a
            a = mid
            sign_a = sign_mid
        else
            b = mid
            sign_b = sign_mid
        end

        root = Arb((a, b))
        verbose && @info "iteration $iteration: $root"
    end

    return (a, b)
end
