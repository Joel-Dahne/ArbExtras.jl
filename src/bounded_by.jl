"""
    bounded_by(f, a::Arf, b::Arf, C::Arf; degree, abs_value, log_bisection, depth_start, maxevals, depth, threaded, verbose)

Return `true` if the function `f` can be shown to be bounded by `C` on
the interval `[a, b]`, i.e. `f(x) <= C` for all `x ∈ [a, b]`,
otherwise return `false`.

This function is similar to first computing the maximum with
[`maximum_enclosure`](@ref) and then check if the computed maximum
satisfies the bound. However if the only thing needed is to check if
the bound holds this method has a number of benefits
- it aborts early if the bound is shown to not hold
- it doesn't try to compute an accurate enclosure of the maximum, it
  only bisects as much as is needed for getting the bound.
- the implementation is simpler and easier to check for correctness.

The maximum of `f` is enclosed by the use of Taylor series through
[`maximum_series`](@ref). The degree of the expansion is kept constant
and the interval is bisected until we can either conclude that the
bound holds on the whole interval or there is subinterval where it
doesn't hold.

The function `f` should support evaluation on both `Arb` and
`ArbSeries` and should return an enclosure of the result in both
cases. The degree used for `ArbSeries` can be set with the `degree`
argument (defaults to 8). If `degree` is negative then it will fall
back to direct evaluation with `Arb` and not make use of
[`maximum_series`](@ref), this is usually much slower but does not
require the function to be implemented for `ArbSeries`.

If `abs_value = true` the instead consider the function `abs(f(x))` on
the interval `[a, b]`.

If `log_bisection = true` then the intervals are bisected in a
logarithmic scale, see [`bisect_interval`](@ref) for details.

The argument `depth_start` bisect the interval using
[`bisect_interval_recursive`](@ref) before starting to compute the
extrema. This can be useful if it is known beforehand that a certain
number of bisections will be necessary before the enclosures get good
enough. It defaults to `0` which corresponds to not bisecting the
interval at all before starting.

The arguments `maxevals` and `depth` can be used to limit the number
of function evaluations and the number of bisections of the interval
respectively. Notice that `depth` takes `depth_start` into account, so
the maximum number of iterations is `depth - depth_start`.

If `threaded = true` then evaluate `f` in parallel on the intervals
using [`Threads.@threads`](@ref).

If `verbose = true` then output information about the process.
"""
function bounded_by(
    f,
    a::Arf,
    b::Arf,
    C::Arf;
    degree::Integer = 8,
    abs_value = false,
    log_bisection = false,
    depth_start::Integer = 0,
    maxevals::Integer = 1000,
    depth::Integer = 20,
    threaded = false,
    verbose = false,
)
    check_interval(a, b)

    maybe_abs = abs_value ? abs : identity

    if a == b
        return maybe_abs(f(Arb(a))) <= C
    end

    # List of intervals
    intervals = bisect_interval_recursive(a, b, depth_start, log_midpoint = log_bisection)

    iterations = 0
    evals = 0

    verbose && @info "iteration: $(lpad(iterations, 2)), " *
          " starting intervals: $(lpad(length(intervals), 3)), "

    while true
        iterations += 1
        evals += length(intervals)

        # Evaluate f on the remaining intervals
        values = similar(intervals, Arb)
        if degree < 0
            if threaded
                Threads.@threads for i in eachindex(intervals)
                    values[i] = maybe_abs(f(Arb(intervals[i])))
                end
            else
                for i in eachindex(intervals)
                    values[i] = maybe_abs(f(Arb(intervals[i])))
                end
            end
        else
            if threaded
                Threads.@threads for i in eachindex(intervals)
                    values[i], _ = maximum_series(f, intervals[i]...; degree, abs_value)
                end
            else
                for (i, (a, b)) in enumerate(intervals)
                    values[i], _ = maximum_series(f, a, b; degree, abs_value)
                end
            end
        end

        # Compute current bounds of maximum for verbose output
        max_low = Arf(-Inf)
        max_upp = Arf(-Inf)

        # Check for each interval if the bound is satisfied
        to_split = falses(length(intervals))
        for i in eachindex(intervals)
            if values[i] > C
                verbose && @info "bound doesn't hold on the interval x = $(intervals[i])"
                verbose && @info "got the maximum of f(x) to be $(values[i])"
                return false
            elseif !(values[i] <= C)
                to_split[i] = true
                if verbose
                    low, upp = getinterval(values[i])
                    Arblib.max!(max_low, max_low, low)
                    Arblib.max!(max_upp, max_upp, upp)
                end
            end
        end
        intervals = bisect_intervals(intervals, to_split, log_midpoint = log_bisection)

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(length(intervals) ÷ 2, 3)), " *
              ifelse(
                  isempty(intervals),
                  "",
                  "maximum on remaining: $(format_interval(max_low, max_upp))",
              )

        non_finite = count(!isfinite, values)
        verbose && non_finite > 0 && @info "non-finite intervals: $non_finite"

        isempty(intervals) && break

        # Check if we have done the maximum number of function
        # evaluations or reached the maximum depth
        if evals >= maxevals || iterations >= depth - depth_start
            if verbose
                evals >= maxevals &&
                    @warn "reached maximum number of evaluations $evals >= $maxevals"
                iterations >= depth - depth_start && @warn "reached maximum depth $depth"
            end
            return false
        end
    end

    return true
end
