"""
    _current_lower_upper_bound(minormax, low, upp, values_low, values_upp, values)

Returns a lower and upper bound of the extremum given that `low` and
`upp` and lower and upper bounds of the extremum and so are all values
in `values_low` and `Values_upp` respectively. It determines
**either** the minimum or the maximum depending on if `minormax = min`
or `minormax = max`.

This is mainly an internal function which implements behaviour which
is common for [`extrema_enclosure`](@ref), [`minimum_enclosure`](@ref)
and [`maximum_enclosure`](@ref).

"""
function _current_lower_upper_bound(
    minormax::Union{typeof(min),typeof(max)},
    low::Arf,
    upp::Arf,
    values_low::Vector{Arf},
    values_upp::Vector{Arf},
    values::Vector{Arb},
)
    if minormax isa typeof(min)
        minormax! = Arblib.min!
    elseif minormax isa typeof(max)
        minormax! = Arblib.max!
    end

    if all(isfinite, values)
        current_low = copy(low)
        current_upp = copy(upp)
        for value in values_low
            minormax!(current_low, current_low, value)
        end
        for value in values_upp
            minormax!(current_upp, current_upp, value)
        end
    elseif minormax isa typeof(min)
        current_low = Arf(-Inf, prec = precision(low))
        current_upp = upp
        for value in values_upp
            isfinite(value) && minormax!(current_upp, current_upp, value)
        end
    elseif minormax isa typeof(max)
        current_low = low
        current_upp = Arf(Inf, prec = precision(upp))
        for value in values_low
            isfinite(value) && minormax!(current_low, current_low, value)
        end
    end

    return current_low, current_upp
end

"""
    extrema_enclosure(f, a::Arf, b::Arf; degree, atol, rtol, lbound_tol, ubound_tol, abs_value, log_bisection, point_value_min, point_value_max, depth_start, maxevals, depth, threaded, verbose)

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

The `atol` and `rtol` argument are used to set the tolerance for when
a subinterval is considered done and not bisected further, see
[`check_tolerance`](@ref) for details. The `lbound_tol` and
`ubound_tol` arguments can optionally be set to give another criteria
for when to stop bisection of subintervals. For the minimum a
subinterval is considered done if the enclosure of its minimum is
strictly greater than `lbound_tol`, for the maximum if its maximum is
strictly less than `ubound_tol`. This can be useful if you mostly care
about whether the minimum/maximum is less/greater than some specific
value, then you don't need to bisect if the current enclosure is good
enough. See also [`bounded_by`](@ref) if you only need to prove that
`f` is bounded by some value.

If `abs_value = true` the compute the extrema of `abs(f(x))` on the
interval `[a, b]`. For the computation of the maximum this is mostly
the same, only difference is that we have to take the absolute value
of the evaluations. For the minimum we have to take into account that
the polynomial might cross zero, in which case the minimum is zero.

If `log_bisection = true` then the intervals are bisected in a
logarithmic scale, see [`bisect_interval`](@ref) for details.

The arguments `point_value_min` and `point_value_max` can optionally
be set to a a priori upper or lower bounds of the min and max
respectively. Typically this is computed by evaluating `f` on some
(possibly well chosen) points before calling this method and taking
the minimum and maximum of the evaluations respectively. This can
allow the method to quicker discard subintervals where the extrema
could not possibly be located.

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
    lbound_tol = Arb(Inf),
    ubound_tol = Arb(-Inf),
    abs_value = false,
    log_bisection = false,
    point_value_min::Arb = Arb(Inf, prec = precision(a)),
    point_value_max::Arb = Arb(-Inf, prec = precision(a)),
    depth_start::Integer = 0,
    maxevals::Integer = 1000,
    depth::Integer = 20,
    threaded = false,
    verbose = false,
)
    check_interval(a, b)

    if a == b
        res = f(Arb(a))
        abs_value && Arblib.abs!(res, res)
        abs_value && Arblib.nonnegative_part!(res, res)
        return res, res
    end

    lbound_tol = convert(Arb, lbound_tol)
    ubound_tol = convert(Arb, ubound_tol)

    # List of intervals
    intervals = bisect_interval_recursive(a, b, depth_start, log_midpoint = log_bisection)

    # Stores the lower and upper bound of the minimum and the maximum
    # on the parts of the interval which are completed. Initially they
    # are set to the neutral element for minimum.
    min_low, min_upp = Arf(Inf, prec = precision(a)), Arf(Inf, prec = precision(a))
    max_low, max_upp = Arf(-Inf, prec = precision(a)), Arf(-Inf, prec = precision(a))

    iterations = 0
    evals = 0

    verbose && @info "iteration: $(lpad(iterations, 2)), " *
          " starting intervals: $(lpad(length(intervals), 3)), "

    while true
        iterations += 1
        evals += length(intervals)

        # Compute enclosure of extrema on each remaining interval
        values_min = similar(intervals, Arb)
        values_max = similar(intervals, Arb)
        if degree < 0
            if threaded
                Threads.@threads for i in eachindex(intervals)
                    v = f(Arb(intervals[i]))
                    abs_value && Arblib.abs!(v, v)
                    values_min[i] = values_max[i] = v
                end
            else
                for i in eachindex(intervals)
                    v = f(Arb(intervals[i]))
                    abs_value && Arblib.abs!(v, v)
                    values_min[i] = values_max[i] = v
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
            for v in point_values
                if isfinite(v)
                    Arblib.min!(point_value_min, point_value_min, v)
                    Arblib.max!(point_value_max, point_value_max, v)
                end
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

        # Compute current lower and upper bound of extrema
        min_current_low, min_current_upp = _current_lower_upper_bound(
            min,
            min_low,
            min_upp,
            values_min_low,
            values_min_upp,
            values_min,
        )
        max_current_low, max_current_upp = _current_lower_upper_bound(
            max,
            max_low,
            max_upp,
            values_max_low,
            values_max_upp,
            values_max,
        )

        min_current_upp = min(min_current_upp, ubound(point_value_min))
        max_current_low = max(max_current_low, lbound(point_value_max))

        # If we are not done split the intervals where extrema could
        # be located and which do not satisfy the tolerance
        to_split = falses(length(intervals))
        for i in eachindex(intervals)
            # Check if the extrema could be located in the interval
            possible_min = values_min_low[i] <= min_current_upp || !isfinite(values_min[i])
            possible_max = max_current_low <= values_max_upp[i] || !isfinite(values_max[i])
            if possible_min || possible_max
                tol_min =
                    check_tolerance(values_min[i]; atol, rtol) ||
                    isfinite(values_min[i]) && values_min[i] > lbound_tol
                tol_max =
                    check_tolerance(values_max[i]; atol, rtol) ||
                    isfinite(values_max[i]) && values_max[i] < ubound_tol

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
                    to_split[i] = true
                end
            end
        end
        intervals = bisect_intervals(intervals, to_split, log_midpoint = log_bisection)

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(length(intervals) ÷ 2, 3)), " *
              "min: $(format_interval(min_current_low, min_current_upp)) " *
              "max: $(format_interval(max_current_low, max_current_upp))"

        if verbose
            non_finite_min = count(!isfinite, values_min)
            non_finite_max = count(!isfinite, values_max)
            (non_finite_min > 0 || non_finite_max > 0) &&
                @info "non-finite intervals: min: $(lpad(non_finite_min, 3)) " *
                      "max: $(lpad(non_finite_max, 3))"
        end

        isempty(intervals) && break

        # Check if we have done the maximum number of function
        # evaluations or reached the maximum depth
        if evals >= maxevals || iterations >= depth - depth_start
            if verbose
                evals >= maxevals &&
                    @info "reached maximum number of evaluations $evals >= $maxevals"
                iterations >= depth - depth_start && @info "reached maximum depth $depth"
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
        Arblib.nonnegative_part!(res_max, res_max)
    end
    return res_min, res_max
end

"""
    minimum_enclosure(f, a::Arf, b::Arf; degree, atol, rtol, lbound_tol, abs_value, log_bisection, point_value_min, depth_start, maxevals, depth, threaded, verbose)

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
    lbound_tol = Arb(Inf),
    abs_value = false,
    log_bisection = false,
    point_value_min::Arb = Arb(Inf, prec = precision(a)),
    depth_start::Integer = 0,
    maxevals::Integer = 1000,
    depth::Integer = 20,
    threaded = false,
    verbose = false,
)
    check_interval(a, b)

    if a == b
        res = f(Arb(a))
        abs_value && Arblib.abs!(res, res)
        abs_value && Arblib.nonnegative_part!(res, res)
        return res
    end

    lbound_tol = convert(Arb, lbound_tol)

    # List of intervals
    intervals = bisect_interval_recursive(a, b, depth_start, log_midpoint = log_bisection)

    # Stores the lower and upper bound of the minimum on the parts of
    # the interval which are completed. Initially they are set to the
    # neutral element for minimum.
    min_low, min_upp = Arf(Inf, prec = precision(a)), Arf(Inf, prec = precision(a))

    iterations = 0
    evals = 0

    verbose && @info "iteration: $(lpad(iterations, 2)), " *
          " starting intervals: $(lpad(length(intervals), 3)), "

    while true
        iterations += 1
        evals += length(intervals)

        # Compute enclosure of minimum on each remaining interval
        values = similar(intervals, Arb)
        if degree < 0
            if threaded
                Threads.@threads for i in eachindex(intervals)
                    v = f(Arb(intervals[i]))
                    abs_value && Arblib.abs!(v, v)
                    values[i] = v
                end
            else
                for i in eachindex(intervals)
                    v = f(Arb(intervals[i]))
                    abs_value && Arblib.abs!(v, v)
                    values[i] = v
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
            for v in point_values
                isfinite(v) && Arblib.min!(point_value_min, point_value_min, v)
            end
        end

        # Get upper and lower bounds of minimum on each interval
        values_low = similar(values, Arf)
        values_upp = similar(values, Arf)
        for i in eachindex(values)
            values_low[i], values_upp[i] = getinterval(values[i])
        end

        # Compute current lower and upper bound of minimum
        min_current_low, min_current_upp = _current_lower_upper_bound(
            min,
            min_low,
            min_upp,
            values_low,
            values_upp,
            values,
        )

        min_current_upp = min(min_current_upp, ubound(point_value_min))

        # If we are not done split the intervals where minimum could
        # be located and which do not satisfy the tolerance
        to_split = falses(length(intervals))
        for i in eachindex(intervals)
            # Check if the minimum could be located in the interval
            if values_low[i] <= min_current_upp || !isfinite(values[i])
                if check_tolerance(values[i]; atol, rtol) ||
                   isfinite(values[i]) && values[i] > lbound_tol
                    # If the interval satisfies the tolerance then add
                    # it to the lower and upper bound of the minimum
                    # for the finished parts of the interval.
                    Arblib.min!(min_low, min_low, values_low[i])
                    Arblib.min!(min_upp, min_upp, values_upp[i])
                else
                    # Otherwise split the interval further
                    to_split[i] = true
                end
            end
        end
        intervals = bisect_intervals(intervals, to_split, log_midpoint = log_bisection)

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(length(intervals) ÷ 2, 3)), " *
              "minimum: $(format_interval(min_current_low, min_current_upp))"

        if verbose
            non_finite = count(!isfinite, values)
            non_finite > 0 && @info "non-finite intervals: $non_finite"
        end

        isempty(intervals) && break

        # Check if we have done the maximum number of function
        # evaluations or reached the maximum depth
        if evals >= maxevals || iterations >= depth - depth_start
            if verbose
                evals >= maxevals &&
                    @info "reached maximum number of evaluations $evals >= $maxevals"
                iterations >= depth - depth_start && @info "reached maximum depth $depth"
            end
            min_low, min_upp = min_current_low, min_current_upp
            break
        end
    end

    res = Arb((min_low, min_upp))
    abs_value && Arblib.nonnegative_part!(res, res)
    return res
end

"""
    maximum_enclosure(f, a::Arf, b::Arf; degree, atol, rtol, ubound_tol, abs_value, log_bisection, point_value_max, depth_start, maxevals, depth, threaded, verbose)

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
    ubound_tol = Arb(-Inf),
    abs_value = false,
    log_bisection = false,
    point_value_max::Arb = Arb(-Inf, prec = precision(a)),
    depth_start::Integer = 0,
    maxevals::Integer = 1000,
    depth::Integer = 20,
    threaded = false,
    verbose = false,
)
    check_interval(a, b)

    if a == b
        res = f(Arb(a))
        abs_value && Arblib.abs!(res, res)
        return res
    end

    ubound_tol = convert(Arb, ubound_tol)

    # List of intervals
    intervals = bisect_interval_recursive(a, b, depth_start, log_midpoint = log_bisection)

    # Stores the lower and upper bound of the maximum on the parts of
    # the interval which are completed. Initially they are set to the
    # neutral element for maximum.
    max_low, max_upp = Arf(-Inf, prec = precision(a)), Arf(-Inf, prec = precision(a))

    iterations = 0
    evals = 0

    verbose && @info "iteration: $(lpad(iterations, 2)), " *
          " starting intervals: $(lpad(length(intervals), 3)), "

    while true
        iterations += 1
        evals += length(intervals)

        # Compute enclosure of maximum on each remaining interval
        values = similar(intervals, Arb)
        if degree < 0
            if threaded
                Threads.@threads for i in eachindex(intervals)
                    v = f(Arb(intervals[i]))
                    abs_value && Arblib.abs!(v, v)
                    values[i] = v
                end
            else
                for i in eachindex(intervals)
                    v = f(Arb(intervals[i]))
                    abs_value && Arblib.abs!(v, v)
                    values[i] = v
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
            for v in point_values
                isfinite(v) && Arblib.max!(point_value_max, point_value_max, v)
            end
        end

        # Get upper and lower bounds of maximum on each interval
        values_low = similar(values, Arf)
        values_upp = similar(values, Arf)
        for i in eachindex(values)
            values_low[i], values_upp[i] = getinterval(Arf, values[i])
        end

        # Compute current lower and upper bound of maximum
        max_current_low, max_current_upp = _current_lower_upper_bound(
            max,
            max_low,
            max_upp,
            values_low,
            values_upp,
            values,
        )

        max_current_low = max(max_current_low, lbound(point_value_max))

        # If we are not done split the intervals where maximum could
        # be located and which do not satisfy the tolerance
        to_split = falses(length(intervals))
        for i in eachindex(intervals)
            # Check if the maximum could be located in the interval
            if max_current_low <= values_upp[i] || !isfinite(values[i])
                if check_tolerance(values[i]; atol, rtol) ||
                   isfinite(values[i]) && values[i] < ubound_tol
                    # If the interval satisfies the tolerance then add
                    # it to the lower and upper bound of the maximum
                    # for the finished parts of the interval.
                    Arblib.max!(max_low, max_low, values_low[i])
                    Arblib.max!(max_upp, max_upp, values_upp[i])
                else
                    # Otherwise split the interval further
                    to_split[i] = true
                end
            end
        end
        intervals = bisect_intervals(intervals, to_split, log_midpoint = log_bisection)

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(length(intervals) ÷ 2, 3)), " *
              "maximum: $(format_interval(max_current_low, max_current_upp))"

        if verbose
            non_finite = count(!isfinite, values)
            non_finite > 0 && @info "non-finite intervals: $non_finite"
        end

        isempty(intervals) && break

        # Check if we have done the maximum number of function
        # evaluations or reached the maximum depth
        if evals >= maxevals || iterations >= depth - depth_start
            if verbose
                evals >= maxevals &&
                    @info "reached maximum number of evaluations $evals >= $maxevals"
                iterations >= depth - depth_start && @info "reached maximum depth $depth"
            end
            max_low, max_upp = max_current_low, max_current_upp
            break
        end
    end

    res = Arb((max_low, max_upp))
    abs_value && Arblib.nonnegative_part!(res, res)
    return res
end
