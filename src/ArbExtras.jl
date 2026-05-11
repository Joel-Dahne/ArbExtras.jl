module ArbExtras

using Arblib, SpecialFunctions

include("compat.jl")

include("temporary.jl")
include("utilities.jl")
include("BisectionLogging.jl")

include("isolate_roots.jl")
include("refine_root.jl")

include("extrema_polynomial.jl")
include("extrema_series.jl")
include("extrema_enclosure.jl")

include("bounded_by.jl")

include("special_functions.jl")

include("integrate.jl")

# Mark primary functions as public (only supported for Julia v1.11 and later)

# Enclosing roots
# isolate_roots.jl
@public isolate_roots
# refine_roots.jl
@public refine_root, refine_root_bisection

# Enclosing extrema
# extrema_polynomial.jl
@public extrema_polynomial, maximum_polynomial, minimum_polynomial
# extrema_series.jl
@public extrema_series, maximum_series, minimum_series, enclosure_series
# extrema_enclosure.jl
@public extrema_enclosure, maximum_enclosure, minimum_enclosure
# bounded_by.jl
@public bounded_by

# Other functionality
# integrate.jl
@public integrate
# utilities.jl
@public bisect_interval, bisect_interval_recursive, bisect_intervals
@public check_tolerance
@public check_interval
@public format_interval
@public taylor_remainder
@public enclosure_ubound, enclosure_lbound, enclosure_getinterval
@public derivative_function
# temporary.jl
@public iscpx
@public compose_zero, compose_zero!

end # module
