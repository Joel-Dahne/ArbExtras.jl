using Arblib, ArbExtras, Test

include("demo_problems.jl")

@testset "ArbExtras" begin
    include("isolate_roots.jl")
    include("refine_root.jl")

    include("extrema_polynomial.jl")
    include("extrema_series.jl")
    include("extrema_enclosure.jl")
end
