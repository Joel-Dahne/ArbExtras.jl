module ArbExtras

using Arblib, SpecialFunctions

include("temporary.jl")
include("utilities.jl")

include("isolate_roots.jl")
include("refine_root.jl")

include("extrema_polynomial.jl")
include("extrema_series.jl")
include("extrema_enclosure.jl")

include("bounded_by.jl")

include("special_functions.jl")

include("integrate.jl")

end # module
