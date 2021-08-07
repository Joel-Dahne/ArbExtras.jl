using Arblib, ArbExtras, Test

@testset "ArbExtras" begin
    include("isolate_roots.jl")
    include("refine_root.jl")

    include("extrema_polynomial.jl")
end
