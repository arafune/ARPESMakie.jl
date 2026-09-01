using ARPESMakie
using Test
using Coverage

@testset "ARPESMakie.jl" begin
    @testset "test for crosshair_heatmap" begin    # Write your tests here.
        include("./crosshair_heatmap.jl")
        include("./waterfall.jl")
        include("./stitch.jl")
    end
end
