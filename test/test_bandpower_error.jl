using CMBLensing
using Measurements
using Test

@testset "binned bandpower standard error" begin
    proj = ProjLambert(Ny=8, Nx=8, θpix=30, T=Float64)
    kwargs = (; ℓedges=[0.0, 1000.0], Cℓfid=ℓ -> 1, err_estimate=true)

    constant_field = FlatFourier(fill(1.0 + 0im, 5, 8), proj)
    constant_error = Measurements.uncertainty.(get_Cℓ(constant_field; kwargs...).Cℓ)
    @test constant_error == [0.0]

    varying_field = FlatFourier(reshape(ComplexF64.(1:40), 5, 8), proj)
    varying_error = Measurements.uncertainty.(get_Cℓ(varying_field; kwargs...).Cℓ)
    @test all(isfinite, varying_error)
    @test all(>=(0), varying_error)
    @test all(>(0), varying_error)
end
