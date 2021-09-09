module ArbExtras

using Arblib

include("temporary.jl")
include("utilities.jl")

include("isolate_roots.jl")
include("refine_root.jl")

include("extrema_polynomial.jl")
include("extrema_series.jl")
include("extrema_enclosure.jl")

end # module
