using Test
using DimensionalData
using DimensionalData: @dim, Dim
using CairoMakie
using ARPES
using ARPES: ARPESData, kx, ky, kz, phi, psi, eV, delay
using ARPESPlots


@testset "ARPESPlots tarpes Heatmap Tests" begin
    data = rand(Float64, 40, 60, 30)
    data[3, 3, 3] = NaN
    data[1, 1, 1] = NaN
    A = ARPESData(
        data,
        (phi(range(-10, 10.0, 40)), eV(range(0, 5.0, 60)), delay(range(-3.0, 5.1, 30))),
    )

    # Heatmap figure generation should not error and return a Figure
    fig = tarpes_evolution_heatmaps(A, 0.0, 0.0)
    @test isa(fig, Figure)
end
