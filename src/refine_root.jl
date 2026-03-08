"""
    refine_root(f, root::Arb; df, atol, rtol, min_iterations, max_iterations, sig_imp_proved, sig_imp_unproved, strict, verbose)
    refine_root(f::ArbPoly, root::Arb; ...)

Refine the given root of the function `f`.

Given a ball `root` containing a unique root of the function `f`, this
method refines the enclosure of the root using successive interval
Newton iterations. The method can only handle simple roots contained
in the interior of the enclosure, and requires the derivative to be
non-zero on the entire enclosure.

The function `f` can be either an `ArbPoly` or a regular function. By
default, the derivative is computed by differentiating the polynomial
in the case of `ArbPoly`, or by using `ArbSeries` in the case of a
regular function. Alternatively, the keyword argument `df` can be set
to a custom function for computing the derivative.

At each iteration, it checks if the required tolerance is met according
to [`check_tolerance`](@ref). The default tolerances are `atol = 0`
and `rtol = 4eps(one(root))`.

If the absolute error does not improve significantly between two
iterations, it stops early, since further iterations are unlikely to
yield substantial improvements. This can happen, for example, when the
function is computed with too low precision. What constitutes a
significant improvement is determined by `sig_imp_proved` if the root
has been proved to exist, and `sig_imp_unproved` otherwise. The
default values are `sig_imp_proved = Mag(1.25)` and `sig_imp_unproved
= Mag(1)`. This check is only performed if at least `min_iterations`
have been completed; this avoids stopping early due to slow
convergence in the beginning. To bypass this check entirely, set
`min_iterations` equal to `max_iterations` (the default is `1`).

If `strict = true`, it will only return an enclosure that has been
proven to contain a unique root within the original enclosure; if it
cannot prove this, it returns `NaN`. If `strict = false`, it will
return an enclosure in all cases. This enclosure is only valid if
`root` indeed contains a unique root, and the returned enclosure may
simply be the starting enclosure `root` itself. This is useful if you
already know that `root` contains a unique root and don't want to
worry about the proof step failing.

If any iteration results in `NaN`, the method returns the enclosure
from the previous iteration (which could be the starting enclosure).
However, if `strict = true` and the existence of a root has not been
proven, it will return `NaN` instead.

If `verbose = true`, it prints the enclosure at each iteration and
provides additional summary information at the end.

In rare cases where `strict = false` and the original enclosure does
not actually contain a root, the method could converge to a potential
root just outside the original enclosure. This happens because
computing the intersection of two balls typically results in an
enclosure slightly larger than the true intersection. This is not an
issue if `strict = true` or if the original enclosure does indeed
contain a root.
"""
function refine_root(
    f,
    root::Arb;
    df = derivative_function(f),
    atol = 0,
    rtol = 4eps(one(root)),
    min_iterations = 1,
    max_iterations = 20,
    sig_imp_proved::Mag = Mag(1.25),
    sig_imp_unproved::Mag = one(Mag),
    strict = true,
    verbose = false,
)
    original_root = root
    root = copy(root)
    mid = midpoint(Arb, root)

    verbose && @info "Starting enclosure: $root"

    err_previous = radius(root)
    isproved = false
    for i = 1:max_iterations
        # Compute new enclosure
        y = f(mid)
        dy = df(root)
        new_root = mid - y / dy

        if !isfinite(new_root)
            verbose && !isfinite(y) && @warn "Non-finite enclosure of y"
            verbose && !isfinite(dy) && @warn "Non-finite enclosure of dy"
            break
        end

        if !Arblib.overlaps(root, new_root)
            if isproved
                error("new enclosure doesn't overlap proved root - f must be wrong.")
            else
                verbose && @warn "New enclosure doesn't overlap root"
                break
            end
        end

        # Note that since Arblib.intersection! only returns an
        # enclosure of the result it's not enough to check that the
        # new enclosure is contained in the interior of previous
        # enclosure. It could happen that the previous enclosure
        # contains numbers not contained in the original enclosure and
        # in that case we are not guaranteed that the root is
        # contained in the original enclosure, only in the previous
        # one.
        if !isproved && Arblib.contains_interior(original_root, new_root)
            verbose && @info "Proved root"
            isproved = true
        end

        Arblib.intersection!(root, root, new_root)
        verbose && @info "Enclosure: $root"

        # If the result satisfies the required tolerance - break
        if isproved && check_tolerance(root; atol, rtol)
            verbose && @info "Tolerance satisfied"
            break
        end

        # If the result did not improve meaningfully compared to the
        # last iteration and we have performed the minimum number of
        # iterations - break. If the root is not yet proved we are
        # more restrictive with breaking early.
        err = radius(root)
        if i >= min_iterations
            err_scaled = isproved ? sig_imp_proved * err : sig_imp_unproved * err

            if err_scaled >= err_previous
                verbose &&
                    @info "Radius only improved from $err_previous to $err - stopping early"
                break
            end
        end
        err_previous = err

        # Update midpoint
        Arblib.set!(mid, Arblib.midref(root))
    end

    if strict && !isproved
        verbose && @warn "could not prove root"
        return Arblib.indeterminate!(root)
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
