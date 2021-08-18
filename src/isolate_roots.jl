"""
    check_interval(f, a::Arf, b::Arf; check_unique = true)

Check if the function `f` has a zero on the interval `[a, b]`.

Returns `(maybe, unique)`. If `maybe` is `false` then `f` is proved to
not have a zero on the interval, if it's `true` then `f` might or
might not have a zero on the interval. If `unique` is `true` then `f`
is guaranteed to have exactly one zero on the interval, if it's
`false` then no information is given.

If `check_unique = false` then don't check for uniqueness, `unique`
will always be `false` in this case.

The function `f` should satisfy the same properties as for
[`isolate_roots`](@ref).

FUTURE WORK:
- Take values on endpoints as arguments to avoid having to compute
  them multiple times.

"""
function check_interval(f, a::Arf, b::Arf; check_unique::Bool = true)
    x = Arb((a, b))
    maybe = Arblib.contains_zero(f(x))

    maybe || return false, false

    check_unique || return true, false

    a_sign = Arblib.sgn_nonzero(f(Arb(a)))
    b_sign = Arblib.sgn_nonzero(f(Arb(b)))

    if a_sign * b_sign < 0
        df = f(ArbSeries([x, 1]))[1]
        return true, !Arblib.contains_zero(df)
    else
        return true, false
    end
end

function check_interval(p::ArbPoly, a::Arf, b::Arf; check_unique::Bool = true)
    x = Arb((a, b))

    check_unique || return Arblib.contains_zero(p(x)), false

    px, dpx = Arblib.evaluate2(p, x)

    Arblib.contains_zero(px) || return false, false

    a_sign = Arblib.sgn_nonzero(p(Arb(a)))
    b_sign = Arblib.sgn_nonzero(p(Arb(b)))

    if a_sign * b_sign < 0
        unique = !Arblib.contains_zero(dpx)
        return true, unique
    elseif a_sign * b_sign > 0
        maybe = Arblib.contains_zero(dpx)
        return maybe, false
    else
        return true, false
    end
end

"""
    isolate_roots(f, a::Arf, b::Arf; depth = 10, check_unique = true, verbose = false)

Isolate the roots of `f` on the interval `[a, b]`.

Returns a tuple `(found, flags)` where `found` is a vector of
subintervals of `[a, b]` given by pairs of `Arf` and `flags` is a
vector of booleans. The output has the following properties:
- The function has no roots on interval outside of the subintervals in
  `found`.
- Subintervals are sorted in increasing order (with no overlap except
  possibly starting and ending with the same point).
- Subintervals with a flag `true` has exactly one (simple) root.
- Subintervals with a flag `false` may or may not contain roots.

If all flags are true then all roots of the function on interval have
been isolated. If there are output subintervals on which the existence
or nonexistence of roots could not be determined, the user may attempt
further searches on those subintervals (possibly with increased
precision and/or increased bounds for the breaking criteria). Note
that roots of multiplicity higher than one and roots located exactly
at endpoints cannot be isolated by the algorithm.

This method closely mimics the function `arb_calc_isolate_roots` in
Arb, both in behaviour and in implementation.

The function `f` should support evaluation on both `Arb` and
`ArbSeries` and should return an enclosure of the result in both
cases. The evaluation on `ArbSeries` is done to compute the derivative
to check uniqueness, this can be skipped by setting `check_unique =
false`, see below. Alternatively `f` can be an `ArbPoly`, in which
case some optimizations are done.

If `check_unique = false` then don't check if roots are unique. In
this case `flags` will always be false and it will simply continue to
split the subintervals until it reaches the maximum depth. The only
benefit to this is that the derivative of `f` doesn't have to be
computed so it can be used when the function of interest doesn't
support evaluating derivatives.

If `verbose = true` then output information after each iteration on
the number of found roots and the number of remaining intervals that
needs to be split further.

FUTURE WORK:
- Parallelize.
- Compute values on endpoints of subintervals once to not have to
  compute them multiple times.
- Should it use `BitVector` or `Vector{Bool}` for flags?
- The current implementation does a breadth-first-search. Consider
  rewriting it to do a depth-first-search instead to reduce memory
  usage in case of many splittings. This might however make it harder
  to parallelize.
"""
function isolate_roots(
    f,
    a::Arf,
    b::Arf;
    depth::Integer = 10,
    check_unique::Bool = true,
    verbose = false,
)
    check_interval(a, b)

    if a == b
        fa = f(Arb(a))
        if Arblib.contains_zero(fa)
            return [(a, b)], [iszero(fa)]
        else
            return Vector{NTuple{2,Arf}}(), Vector{Bool}()
        end
    end

    intervals = [(a, b)]
    iterations = 0

    found = Vector{NTuple{2,Arf}}()
    flags = BitVector()

    while !isempty(intervals) && iterations < depth
        iterations += 1

        next_intervals = Vector{eltype(intervals)}()

        for interval in intervals
            maybe, unique = check_interval(f, interval...; check_unique)

            if unique
                push!(found, interval)
                push!(flags, true)
            elseif maybe
                if iterations < depth
                    push!(next_intervals, bisect_interval(interval...)...)
                else
                    # If we are on the last iteration don't split the interval
                    push!(next_intervals, interval)
                end
            end
        end

        intervals = next_intervals

        verbose &&
            iterations < depth &&
            @info "isolate_roots iteration: $(lpad(iterations, 2)), " *
                  "found: $(lpad(length(found), 2)), " *
                  "remaining intervals: $(length(intervals) ÷ 2)"
    end

    verbose &&
        iterations == depth &&
        @info "isolate_roots iteration: $(lpad(iterations, 2)), " *
              "found: $(lpad(length(found), 2)), " *
              "remaining intervals: $(length(intervals))"

    found = [found; intervals]
    flags = [flags; zeros(Bool, length(intervals))]

    p = sortperm(found, by = interval -> interval[1])

    return found[p], flags[p]
end
