abstract type BireStaticOp{T} <: StaticOpWithAdjoint{T} end

struct BireStatic{T} <: BireStaticOp{T}
    α::Field
end

BireStatic(α::Field) = BireStatic{real(eltype(α))}(α)

function (R::BireStatic)(f::Field)
    α_map = Map(R.α)
    f_iqu = IQUMap(f)

    I = f_iqu.I
    Q = f_iqu.Q
    U = f_iqu.U

    cos2α = cos.(2 .* α_map)
    sin2α = sin.(2 .* α_map)

    Q_rot = cos2α .* Q .- sin2α .* U
    U_rot = sin2α .* Q .+ cos2α .* U

    f_iqu_rot = copy(f_iqu)
    f_iqu_rot.arr[:, :, 1] .= I.arr
    f_iqu_rot.arr[:, :, 2] .= Q_rot.arr
    f_iqu_rot.arr[:, :, 3] .= U_rot.arr

    return IEBFourier(f_iqu_rot)
end

import Base: *, \
import Zygote

*(R::BireStatic, f::Field) = R(f)

# Override the generic StaticOp adjoint so alpha gradients flow through
# birefringence multiplies via the custom ChainRules rrules below.
Zygote.@adjoint function *(R::BireStatic, f::Field)
    y, back = rrule(*, R, f)
    function pullback(Δ)
        _, ∂R, ∂f = back(Δ)
        return (∂R, ∂f)
    end
    return y, pullback
end

# adjoint = rotation by -α
*(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) = BireStatic(-Radj.parent.α)(f)

Zygote.@adjoint function *(Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
    y, back = rrule(*, Radj, f)
    function pullback(Δ)
        _, ∂Radj, ∂f = back(Δ)
        return (∂Radj, ∂f)
    end
    return y, pullback
end

# inverse
\(R::BireStatic, f::Field) = BireStatic(-R.α)(f)

# inverse adjoint
\(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) = BireStatic(Radj.parent.α)(f)

import ChainRulesCore
const CRC = ChainRulesCore
import ChainRulesCore: rrule, NoTangent, ZeroTangent, Tangent, unthunk

# rrule for R * f
function rrule(::typeof(*), R::BireStatic, f::Field)
    α_map = Map(R.α)
    f_iqu = IQUMap(f)

    I = f_iqu.I
    Q = f_iqu.Q
    U = f_iqu.U

    cos2α = cos.(2 .* α_map)
    sin2α = sin.(2 .* α_map)

    Q_rot = cos2α .* Q .- sin2α .* U
    U_rot = sin2α .* Q .+ cos2α .* U

    f_iqu_rot = copy(f_iqu)
    f_iqu_rot.arr[:, :, 1] .= I.arr
    f_iqu_rot.arr[:, :, 2] .= Q_rot.arr
    f_iqu_rot.arr[:, :, 3] .= U_rot.arr
    y = IEBFourier(f_iqu_rot)

    function pullback(Δ)
        Δu = unthunk(Δ)
        Δ_iqu = IQUMap(Δu)

        ΔQ = Δ_iqu.Q
        ΔU = Δ_iqu.U

        # ∂L/∂f = R' Δ = rotation by -α
        Q_back =  cos2α .* ΔQ .+ sin2α .* ΔU
        U_back = -sin2α .* ΔQ .+ cos2α .* ΔU

        f_iqu_back = copy(f_iqu)
        f_iqu_back.arr[:, :, 1] .= zero(I.arr)
        f_iqu_back.arr[:, :, 2] .= Q_back.arr
        f_iqu_back.arr[:, :, 3] .= U_back.arr
        ∂f = IEBFourier(f_iqu_back)

        # dQ_rot/dα = -2 sin(2α) Q - 2 cos(2α) U
        # dU_rot/dα =  2 cos(2α) Q - 2 sin(2α) U
        ∂α_from_Q = -2 .* sin2α .* Q .- 2 .* cos2α .* U
        ∂α_from_U =  2 .* cos2α .* Q .- 2 .* sin2α .* U

        ∂α_map = similar(α_map)
        ∂α_map.arr .= ∂α_from_Q.arr .* ΔQ.arr .+ ∂α_from_U.arr .* ΔU.arr

        ∂R = Tangent{typeof(R)}(α = ∂α_map)
        return (NoTangent(), ∂R, ∂f)
    end

    return y, pullback
end

# rrule for R' * f
function rrule(::typeof(*), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
    R = Radj.parent
    α_map = Map(R.α)
    f_iqu = IQUMap(f)

    I = f_iqu.I
    Q = f_iqu.Q
    U = f_iqu.U

    cos2α = cos.(2 .* α_map)
    sin2α = sin.(2 .* α_map)

    # R' = rotation by -α
    Q_rot =  cos2α .* Q .+ sin2α .* U
    U_rot = -sin2α .* Q .+ cos2α .* U

    f_iqu_rot = copy(f_iqu)
    f_iqu_rot.arr[:, :, 1] .= I.arr
    f_iqu_rot.arr[:, :, 2] .= Q_rot.arr
    f_iqu_rot.arr[:, :, 3] .= U_rot.arr
    y = IEBFourier(f_iqu_rot)

    function pullback(Δ)
        Δu = unthunk(Δ)
        Δ_iqu = IQUMap(Δu)

        ΔQ = Δ_iqu.Q
        ΔU = Δ_iqu.U

        # ∂L/∂f = R Δ = rotation by +α
        Q_back = cos2α .* ΔQ .- sin2α .* ΔU
        U_back = sin2α .* ΔQ .+ cos2α .* ΔU

        f_iqu_back = copy(f_iqu)
        f_iqu_back.arr[:, :, 1] .= zero(I.arr)
        f_iqu_back.arr[:, :, 2] .= Q_back.arr
        f_iqu_back.arr[:, :, 3] .= U_back.arr
        ∂f = IEBFourier(f_iqu_back)

        # derivative wrt α for adjoint operator
        ∂α_from_Q =  2 .* sin2α .* Q .- 2 .* cos2α .* U
        ∂α_from_U =  2 .* cos2α .* Q .+ 2 .* sin2α .* U

        ∂α_map = similar(α_map)
        ∂α_map.arr .= ∂α_from_Q.arr .* ΔQ.arr .+ ∂α_from_U.arr .* ΔU.arr

        ∂Radj = Tangent{typeof(Radj)}(
            parent = Tangent{typeof(R)}(α = ∂α_map)
        )

        return (NoTangent(), ∂Radj, ∂f)
    end

    return y, pullback
end

# constructor rules
function rrule(::Type{<:BireStatic}, α::Field)
    y = BireStatic(α)
    function pullback(Δ)
        Δu = unthunk(Δ)
        ∂α = hasproperty(Δu, :α) ? getproperty(Δu, :α) : ZeroTangent()
        return (NoTangent(), ∂α)
    end
    return y, pullback
end

function rrule(::Type{<:BireStatic{T}}, α::Field) where {T}
    y = BireStatic{T}(α)
    function pullback(Δ)
        Δu = unthunk(Δ)
        ∂α = hasproperty(Δu, :α) ? getproperty(Δu, :α) : ZeroTangent()
        return (NoTangent(), ∂α)
    end
    return y, pullback
end