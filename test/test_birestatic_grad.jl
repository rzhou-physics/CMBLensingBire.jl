# Finite-difference checks for the α-gradient of the birefringence operator
# `BireStatic`, for BOTH  R*f  and  R\f.
#
# The `R\f` path is the one exercised by `unmix` (dataset.jl) and hence by the
# ϕ-α step of `MAP_joint`. Before the `\` adjoints were added, that gradient was
# silently dropped, biasing the joint birefringence estimate. Tests 2-4 fail on
# the old code and pass once the `\` rules are in place.
#
# This precompiles CMBLensing, so run it on a compute node, e.g.
#   julia --project=. test/test_birestatic_grad.jl
# (switch storage=CuArray below if you want to check the GPU path too).

using CMBLensing, PythonPlot, Zygote, LinearAlgebra, Random, Printf
using Distributions: logpdf

# ---- small CPU sim to obtain consistent (f, ϕ, α, ds) -----------------------
sim = load_sim(
    seed     = 1,
    Cℓ       = camb(r=0.001, ℓmax=6000),
    Cℓn      = noiseCℓs(μKarcminT=1.0),
    θpix     = 3,
    T        = Float64,
    Nside    = 16,
    pol      = :IP,
    ℓmin_ϕ   = 100, ℓmax_ϕ = 2000, Nbins_ϕ = 2,
    ℓmin_α   = 100, ℓmax_α = 2000, Nbins_α = 2,
    qe_noise = false,
    storage  = Array,
)
f, ϕ, α, ds = sim.f, sim.ϕ, sim.α, sim.ds

Random.seed!(2)

# Directional finite-difference check of an AD gradient.
# `L`  : scalar function of a spin-0 field
# `x0` : base point (a spin-0 Map/field)
function fd_check(name, L, x0; ε=1e-6, tol=1e-4)
    δ = zero(Map(x0))
    δ.arr .= randn(eltype(δ.arr), size(δ.arr))
    g       = Zygote.gradient(L, x0)[1]              # AD gradient
    ad      = real(dot(g, δ))                        # directional derivative (AD)
    fd      = (L(x0 + ε*δ) - L(x0 - ε*δ)) / (2ε)     # directional derivative (FD)
    relerr  = abs(ad - fd) / max(abs(fd), eps())
    status  = relerr < tol ? "OK" : "*** FAIL ***"
    @printf("%-30s  AD=% .6e  FD=% .6e  relerr=%.2e  %s\n", name, ad, fd, relerr, status)
    return relerr < tol
end

# Fixed random cotangent-carrying weight, in the output (IEB) field type.
# (A norm loss ½‖R f‖² is α-independent because rotations are orthogonal, so we
#  contract against a fixed random field instead.)
wtmpl = BireStatic(α) * f
w = zero(wtmpl); w.arr .= randn(eltype(w.arr), size(w.arr))

ok = Bool[]

# Test 1:  R(α) * f      (was already correct)
push!(ok, fd_check("R(α) * f",          a -> real(dot(w, BireStatic(a) * f)), α))
# Test 2:  R(α) \ f      (the fix — dropped before; this is the unmix path)
push!(ok, fd_check("R(α) \\ f  [unmix]", a -> real(dot(w, BireStatic(a) \ f)), α))
# Test 3:  R(α)' \ f     (adjoint-inverse, for completeness)
push!(ok, fd_check("R(α)' \\ f",         a -> real(dot(w, BireStatic(a)' \ f)), α))

# Test 4 (integration):  α°-gradient of the mixed-parametrization logpdf, which
# internally calls unmix -> R(α)\ . This is exactly the MAP_joint ϕ-α step path.
dsθ = copy(ds())
dsθ.G = I
dsθ.J = I                                    # mirror MAP_joint  => α° == α
m  = CMBLensing.mix(dsθ; f=f, ϕ=ϕ, α=α)
f°, ϕ°, α°0 = m.f°, m.ϕ°, m.α°
Lα°(a°) = logpdf(CMBLensing.Mixed(dsθ); f°=f°, ϕ°=ϕ°, α°=a°, θ=(;))
push!(ok, fd_check("logpdf(Mixed) wrt α°", Lα°, α°0))

println()
if all(ok)
    println("ALL BireStatic gradient checks PASSED ✓")
else
    error("Some BireStatic gradient checks FAILED — see *** FAIL *** rows above.")
end
