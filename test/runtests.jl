using Arblib, ArbExtras, Random, SpecialFunctions, Test

include("demo_problems.jl")

@testset "ArbExtras" begin
    include("utilities.jl")

    include("isolate_roots.jl")
    include("refine_root.jl")

    include("extrema_polynomial.jl")
    include("extrema_series.jl")
    include("extrema_enclosure.jl")

    include("special_functions.jl")
end
