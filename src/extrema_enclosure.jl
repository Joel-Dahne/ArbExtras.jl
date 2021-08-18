"""
    extrema_enclosure(f, a::Arf, b::Arf; degree, atol, rtol, abs_value, point_value_min, point_value_max, maxevals, depth, threaded, verbose)

Compute both the minimum and maximum of the function `f` on the
interval `[a, b]` and return them as a 2-tuple.

The extrema are computed by enclosing the extrema of the function on
an interval using Taylor series through [`extrema_series`](@ref). The
degree of the expansion is kept constant and the interval is bisected
until the required tolerance is met.

The function `f` should support evaluation on both `Arb` and
`ArbSeries` and should return an enclosure of the result in both
cases. The degree used for `ArbSeries` can be set with the `degree`
argument (defaults to 8). If `degree` is negative then it will fall
back to direct evaluation with `Arb` and not make use of
[`extrema_series`](@ref), this is usually much slower but does not
require the function to be implemented for `ArbSeries`.

If `abs_value = true` the compute the extrema of `abs(f(x))` on the
interval `[a, b]`. For the computation of the maximum this is mostly
the same, only difference is that we have to take the absolute value
of the evaluations. For the minimum we have to take into account that
the polynomial might cross zero, in which case the minimum is zero.

The arguments `point_value_min` and `point_value_max` can optionally
be set to a a priori upper or lower bounds of the min and max
respectively. Typically this is computed by evaluating `f` on some
(possibly well chosen) points before calling this method and taking
the minimum and maximum of the evaluations respectively. This can
allow the method to quicker discard subintervals where the extrema
could not possibly be located.

The arguments `maxevals` and `depth` can be used to limit the number
of function evaluations and the number of bisections of the interval
respectively.

If `threaded = true` then evaluate `f` in parallel on the intervals
using [`Threads.@threads`](@ref).

If `verbose = true` then output information about the process.

TODO: Currently this always computes both minimum and maximum on all
subintervals. Change it so that it takes into account if the
minimum/maximum could be located in the subinterval or not.

"""
function extrema_enclosure(
    f,
    a::Arf,
    b::Arf;
    degree::Integer = 8,
    atol = 0,
    rtol = sqrt(eps(one(a))),
    abs_value = false,
    point_value_min::Arb = Arb(Inf, prec = precision(a)),
    point_value_max::Arb = Arb(-Inf, prec = precision(a)),
    maxevals::Integer = 1000,
    depth::Integer = 20,
    threaded = false,
    verbose = false,
)
    check_interval(a, b)

    maybe_abs = abs_value ? abs : identity

    if a == b
        fa = maybe_abs(f(Arb(a)))
        abs_value && Arblib.nonnegative_part!(fa, fa)
        return fa, fa
    end

    # List of intervals left to process
    intervals = [(a, b)]

    # Stores the lower and upper bound of the minimum and the maximum
    # on the parts of the interval which are completed. Initially they
    # are set to the neutral element for minimum.
    min_low, min_upp = Arf(Inf, prec = precision(a)), Arf(Inf, prec = precision(a))
    max_low, max_upp = Arf(-Inf, prec = precision(a)), Arf(-Inf, prec = precision(a))

    iterations = 0
    evals = 0

    while true
        iterations += 1
        evals += length(intervals)

        # Compute enclosure of extrema on each remaining interval
        values_min = similar(intervals, Arb)
        values_max = similar(intervals, Arb)
        if degree < 0
            if threaded
                Threads.@threads for i in eachindex(intervals)
                    values_min[i] = values_max[i] = maybe_abs(f(Arb(intervals[i])))
                end
            else
                for i in eachindex(intervals)
                    values_min[i] = values_max[i] = maybe_abs(f(Arb(intervals[i])))
                end
            end
        else
            point_values = similar(values_min)
            if threaded
                Threads.@threads for i in eachindex(intervals)
                    values_min[i], values_max[i], point_values[i] =
                        extrema_series(f, intervals[i]...; degree, abs_value)
                end
            else
                for (i, (a, b)) in enumerate(intervals)
                    values_min[i], values_max[i], point_values[i] =
                        extrema_series(f, a, b; degree, abs_value)
                end
            end
            for point_value in filter(isfinite, point_values)
                Arblib.min!(point_value_min, point_value_min, point_value)
                Arblib.max!(point_value_max, point_value_max, point_value)
            end
        end

        # Get upper and lower bounds of extrema on each interval
        values_min_low = similar(values_min, Arf)
        values_min_upp = similar(values_min, Arf)
        values_max_low = similar(values_max, Arf)
        values_max_upp = similar(values_max, Arf)
        for i in eachindex(intervals)
            values_min_low[i], values_min_upp[i] = getinterval(values_min[i])
            values_max_low[i], values_max_upp[i] = getinterval(values_max[i])
        end

        # Compute current lower and upper bound of extrema on both the
        # completed parts of the interval and the remaining ones.
        if all(isfinite, values_min)
            min_current_low = min(min_low, minimum(values_min_low))
            min_current_upp = min(min_upp, minimum(values_min_upp))
        else
            min_current_low = Arf(-Inf, prec = precision(a))
            min_current_upp = min(
                min_upp,
                minimum(
                    filter(isfinite, values_min_upp),
                    init = Arf(Inf, prec = precision(a)),
                ),
            )
        end
        if all(isfinite, values_max)
            max_current_low = max(max_low, maximum(values_max_low))
            max_current_upp = max(max_upp, maximum(values_max_upp))
        else
            max_current_low = max(
                max_low,
                maximum(
                    filter(isfinite, values_max_low),
                    init = Arf(-Inf, prec = precision(a)),
                ),
            )
            max_current_upp = Arf(Inf, prec = precision(a))
        end

        min_current_upp = min(min_current_upp, ubound(point_value_min))
        max_current_low = max(max_current_low, lbound(point_value_max))

        # If we are not done split the intervals where extrema could
        # be located and which do not satisfy the tolerance
        next_intervals = Vector{eltype(intervals)}()
        for i in eachindex(intervals)
            # Check if the extrema could be located in the interval
            possible_min = values_min_low[i] <= min_current_upp || !isfinite(values_min[i])
            possible_max = max_current_low <= values_max_upp[i] || !isfinite(values_max[i])
            if possible_min || possible_max
                tol_min = check_tolerance(values_min[i]; atol, rtol)
                tol_max = check_tolerance(values_max[i]; atol, rtol)

                if (!possible_min || tol_min) && (!possible_max || tol_max)
                    # If the interval satisfies the tolerance then add
                    # it to the lower and upper bound of the extrema
                    # for the finished parts of the interval.
                    Arblib.min!(min_low, min_low, values_min_low[i])
                    Arblib.min!(min_upp, min_upp, values_min_upp[i])
                    Arblib.max!(max_low, max_low, values_max_low[i])
                    Arblib.max!(max_upp, max_upp, values_max_upp[i])
                else
                    # Otherwise split the interval further
                    push!(next_intervals, bisect_interval(intervals[i]...)...)
                end
            end
        end
        intervals = next_intervals

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(length(intervals) ÷ 2, 3)), " *
              "minimum: $(Float64.((min_current_low, min_current_upp)))" *
              "maximum: $(Float64.((max_current_low, max_current_upp)))"

        isempty(intervals) && break

        # Check if we have done the maximum number of function
        # evaluations or reached the maximum depth
        if evals >= maxevals || iterations >= depth
            if verbose
                evals >= maxevals &&
                    @info "reached maximum number of evaluations $evals >= $maxevals"
                iterations >= depth && @info "reached maximum depth $depth"
            end
            min_low, min_upp = min_current_low, min_current_upp
            max_low, max_upp = max_current_low, max_current_upp
            break
        end
    end

    res_min = Arb((min_low, min_upp))
    res_max = Arb((max_low, max_upp))
    if abs_value
        Arblib.nonnegative_part!(res_min, res_min)
    end
    return res_min, res_max
end

"""
    minimum_enclosure(f, a::Arf, b::Arf; degree, atol, rtol, abs_value, point_value_min, maxevals, depth, threaded, verbose)

Compute the minimum of the function `f` on the interval `[a, b]`.

Takes the same arguments as [`extrema_polynomial`](@ref). The
algorithm is also the same except that it only looks for the minimum.
"""
function minimum_enclosure(
    f,
    a::Arf,
    b::Arf;
    degree::Integer = 8,
    atol = 0,
    rtol = sqrt(eps(one(a))),
    abs_value = false,
    point_value_min::Arb = Arb(Inf, prec = precision(a)),
    maxevals::Integer = 1000,
    depth::Integer = 20,
    threaded = false,
    verbose = false,
)
    check_interval(a, b)

    maybe_abs = abs_value ? abs : identity

    if a == b
        res = maybe_abs(f(Arb(a)))
        if abs_value
            return Arblib.nonnegative_part!(res, res)
        else
            return res
        end
    end

    # List of intervals left to process
    intervals = [(a, b)]

    # Stores the lower and upper bound of the minimum on the parts of
    # the interval which are completed. Initially they are set to the
    # neutral element for minimum.
    min_low, min_upp = Arf(Inf, prec = precision(a)), Arf(Inf, prec = precision(a))

    iterations = 0
    evals = 0
    while true
        iterations += 1
        evals += length(intervals)

        # Compute enclosure of minimum on each remaining interval
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
            point_values = similar(values)
            if threaded
                Threads.@threads for i in eachindex(intervals)
                    values[i], point_values[i] =
                        minimum_series(f, intervals[i]...; degree, abs_value)
                end
            else
                for (i, (a, b)) in enumerate(intervals)
                    values[i], point_values[i] = minimum_series(f, a, b; degree, abs_value)
                end
            end
            for point_value in filter(isfinite, point_values)
                Arblib.min!(point_value_min, point_value_min, point_value)
            end
        end

        # Get upper and lower bounds of minimum on each interval
        values_low = similar(values, Arf)
        values_upp = similar(values, Arf)
        for i in eachindex(values)
            values_low[i], values_upp[i] = getinterval(values[i])
        end

        # Compute current lower and upper bound of minimum on both the
        # completed parts of the interval and the remaining ones.
        if all(isfinite, values)
            min_current_low = min(min_low, minimum(values_low))
            min_current_upp = min(min_upp, minimum(values_upp))
        else
            min_current_low = Arf(-Inf, prec = precision(a))
            min_current_upp = min(
                min_upp,
                minimum(filter(isfinite, values_upp), init = Arf(Inf, prec = precision(a))),
            )
        end

        min_current_upp = min(min_current_upp, ubound(point_value_min))

        # If we are not done split the intervals where minimum could
        # be located and which do not satisfy the tolerance
        next_intervals = Vector{eltype(intervals)}()
        for i in eachindex(intervals)
            # Check if the minimum could be located in the interval
            if values_low[i] <= min_current_upp || !isfinite(values[i])
                if check_tolerance(values[i]; atol, rtol)
                    # If the interval satisfies the tolerance then add
                    # it to the lower and upper bound of the minimum
                    # for the finished parts of the interval.
                    Arblib.min!(min_low, min_low, values_low[i])
                    Arblib.min!(min_upp, min_upp, values_upp[i])
                else
                    # Otherwise split the interval further
                    push!(next_intervals, bisect_interval(intervals[i]...)...)
                end
            end
        end
        intervals = next_intervals

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(length(intervals) ÷ 2, 3)), " *
              "minimum: $(Float64.((min_current_low, min_current_upp)))"

        isempty(intervals) && break

        # Check if we have done the maximum number of function
        # evaluations or reached the maximum depth
        if evals >= maxevals || iterations >= depth
            if verbose
                evals >= maxevals &&
                    @info "reached maximum number of evaluations $evals >= $maxevals"
                iterations >= depth && @info "reached maximum depth $depth"
            end
            min_low, min_upp = min_current_low, min_current_upp
            break
        end
    end

    res = Arb((min_low, min_upp))
    if abs_value
        Arblib.nonnegative_part!(res, res)
    end
    return res
end

"""
    maximum_enclosure(f, a::Arf, b::Arf; degree, atol, rtol, abs_value, point_value_max, maxevals, depth, threaded, verbose)

Compute the maximum of the function `f` on the interval `[a, b]`.

Takes the same arguments as [`extrema_polynomial`](@ref). The
algorithm is also the same except that it only looks for the maximum.
"""
function maximum_enclosure(
    f,
    a::Arf,
    b::Arf;
    degree::Integer = 8,
    atol = 0,
    rtol = sqrt(eps(one(a))),
    abs_value = false,
    point_value_max::Arb = Arb(-Inf, prec = precision(a)),
    maxevals::Integer = 1000,
    depth::Integer = 20,
    threaded = false,
    verbose = false,
)
    check_interval(a, b)

    maybe_abs = abs_value ? abs : identity

    if a == b
        return maybe_abs(f(Arb(a)))
    end

    # List of intervals left to process
    intervals = [(a, b)]

    # Stores the lower and upper bound of the maximum on the parts of
    # the interval which are completed. Initially they are set to the
    # neutral element for maximum.
    max_low, max_upp = Arf(-Inf, prec = precision(a)), Arf(-Inf, prec = precision(a))

    iterations = 0
    evals = 0
    while true
        iterations += 1
        evals += length(intervals)

        # Compute enclosure of maximum on each remaining interval
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
            point_values = similar(values)
            if threaded
                Threads.@threads for i in eachindex(intervals)
                    values[i], point_values[i] =
                        maximum_series(f, intervals[i]...; degree, abs_value)
                end
            else
                for (i, (a, b)) in enumerate(intervals)
                    values[i], point_values[i] = maximum_series(f, a, b; degree, abs_value)
                end
            end
            for point_value in filter(isfinite, point_values)
                Arblib.max!(point_value_max, point_value_max, point_value)
            end
        end

        # Get upper and lower bounds of maximum on each interval
        values_low = similar(values, Arf)
        values_upp = similar(values, Arf)
        for i in eachindex(values)
            values_low[i], values_upp[i] = getinterval(Arf, values[i])
        end

        # Compute current lower and upper bound of maximum on both the
        # completed parts of the interval and the remaining ones.
        if all(isfinite, values)
            max_current_low = max(max_low, maximum(values_low))
            max_current_upp = max(max_upp, maximum(values_upp))
        else
            max_current_low = max(
                max_low,
                maximum(
                    filter(isfinite, values_low),
                    init = Arf(-Inf, prec = precision(a)),
                ),
            )
            max_current_upp = Arf(Inf, prec = precision(a))
        end

        max_current_low = max(max_current_low, lbound(point_value_max))

        # If we are not done split the intervals where maximum could
        # be located and which do not satisfy the tolerance
        next_intervals = Vector{eltype(intervals)}()
        for i in eachindex(intervals)
            # Check if the maximum could be located in the interval
            if max_current_low <= values_upp[i] || !isfinite(values[i])
                if check_tolerance(values[i]; atol, rtol)
                    # If the interval satisfies the tolerance then add
                    # it to the lower and upper bound of the maximum
                    # for the finished parts of the interval.
                    Arblib.max!(max_low, max_low, values_low[i])
                    Arblib.max!(max_upp, max_upp, values_upp[i])
                else
                    # Otherwise split the interval further
                    push!(next_intervals, bisect_interval(intervals[i]...)...)
                end
            end
        end
        intervals = next_intervals

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(length(intervals) ÷ 2, 3)), " *
              "maximum: $(Float64.((max_current_low, max_current_upp)))"

        isempty(intervals) && break

        # Check if we have done the maximum number of function
        # evaluations or reached the maximum depth
        if evals >= maxevals || iterations >= depth
            if verbose
                evals >= maxevals &&
                    @info "reached maximum number of evaluations $evals >= $maxevals"
                iterations >= depth && @info "reached maximum depth $depth"
            end
            max_low, max_upp = max_current_low, max_current_upp
            break
        end
    end

    return Arb((max_low, max_upp))
end
