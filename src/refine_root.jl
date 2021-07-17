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

At each iteration it checks if the required tolerance is met. The
absolute error is given by the diameter of the enclosure and the
relative error is given by the absolute error divided by the
enclosure. The default tolerances are `atol = 0` and `rtol =
4eps(one(root))`.

If the absolute error improve by at least a factor 1.5 between two
iterations then it stops early since it's unlikely that further
iterations would improve the result much. This can for example happen
when the function is computed with too low precision. This check is
only done if `min_iterations` have been performed, this is to avoid
stopping early due to slow convergence in the beginning. To avoid this
check entirely `min_iterations` can be put to the same as
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
        series = ArbSeries([root, 1])
    end

    verbose && @info "enclosure: $root"

    error_previous = 2Arblib.radius(Arb, root)
    isproved = false
    for i = 1:max_iterations
        # Compute new enclosure
        y = f(mid)
        if f isa ArbPoly
            dy = df(root)
        else
            series
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
        error = 2Arblib.radius(Arb, root)
        if error <= atol || error / abs(root) <= rtol
            verbose && @info "tolerance satisfied"
            break
        end

        # If the result did not improve compared to the last iteration
        # and we have performed the minimum number of iterations -
        # break
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
