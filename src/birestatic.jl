# abstract type BireStaticOp{T} <: StaticOpWithAdjoint{T} end

# struct BireStatic{T} <: BireStaticOp{T}
#     α::Field
# end

# BireStatic(α::Field) = BireStatic{real(eltype(α))}(α)

# # Forward operator R(α)*f
# function (R::BireStatic)(f::Field)
#     α = Map(R.α)
#     f_iqu = IQUMap(f)

#     Q = f_iqu.Q
#     U = f_iqu.U

#     Q_rot =  cos.(2 .* α) .* Q .- sin.(2 .* α) .* U
#     U_rot =  sin.(2 .* α) .* Q .+ cos.(2 .* α) .* U

#     f_iqu_rot = copy(f_iqu)
#     f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
#     f_iqu_rot.arr[:, :, end]   .= U_rot.arr

#     return IQUFourier(f_iqu_rot)
# end

# import Base: *
# *(R::BireStatic, f::Field) = R(f)

# # Adjoint operator R' = R(-α)
# import Base: adjoint
# *(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) =
#     BireStatic(-Radj.parent.α)(f)

# # Inverse operator R \ f = R(-α)*f
# import Base: \
# \(R::BireStatic, f::Field) =
#     BireStatic(-R.α)(f)

# # Inverse-adjoint operator R' \ f = R(+α)*f
# \(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) =
#     BireStatic(Radj.parent.α)(f)

# # AD rule for R*f: no gradients, treat α as constant
# import ChainRulesCore: rrule, NoTangent, ZeroTangent, AbstractTangent
# function rrule(::typeof(*), R::BireStatic, f::Field)
#     y = R * f
#     pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
#     return y, pullback
# end

# # AD rule for R'*f
# function rrule(::typeof(*), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
#     y = Radj * f
#     pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
#     return y, pullback
# end

# # AD rule for R\f
# function rrule(::typeof(\), R::BireStatic, f::Field)
#     y = R \ f
#     pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
#     return y, pullback
# end

# # AD rule for R'\f
# function rrule(::typeof(\), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
#     y = Radj \ f
#     pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
#     return y, pullback
# end

# # Kill any tangent inputs: R*t → 0
# *(R::BireStatic, t::AbstractTangent) = ZeroTangent()
# *(Radj::Adjoint{<:Any,<:BireStatic}, t::AbstractTangent) = ZeroTangent()
# \(R::BireStatic, t::AbstractTangent) = ZeroTangent()
# \(Radj::Adjoint{<:Any,<:BireStatic}, t::AbstractTangent) = ZeroTangent()

abstract type BireStaticOp{T} <: StaticOpWithAdjoint{T} end

struct BireStatic{T} <: BireStaticOp{T}
    α::Field
end

BireStatic(α::Field) = BireStatic{real(eltype(α))}(α)

function (R::BireStatic)(f::Field)
    α = Map(R.α)
    f_iqu = IQUMap(f)

    Q = f_iqu.Q
    U = f_iqu.U

    Q_rot =  cos.(2 .* α) .* Q .- sin.(2 .* α) .* U
    U_rot =  sin.(2 .* α) .* Q .+ cos.(2 .* α) .* U

    f_iqu_rot = copy(f_iqu)
    f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
    f_iqu_rot.arr[:, :, end]   .= U_rot.arr

    return IQUFourier(f_iqu_rot)
end

import Base: *
*(R::BireStatic, f::Field) = R(f)

# -----------------------------------------
#  Adjoint operator:  R' * f = R(-α) * f
# -----------------------------------------
import Base: adjoint

*(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) =
    BireStatic(-Radj.parent.α)(f)


# -----------------------------------------
#  Inverse operator:  R \ f = R(-α) * f
# -----------------------------------------
import Base: \

\(R::BireStatic, f::Field) =
    BireStatic(-R.α)(f)


# ------------------------------------------------------
#  Inverse-adjoint:  R' \ f = R(+α) * f
# ------------------------------------------------------
\(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) =
    BireStatic(Radj.parent.α)(f)

import ChainRulesCore: rrule, NoTangent, ZeroTangent

function rrule(::typeof(*), R::BireStatic, f::Field)
    y = R * f
    pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
    return y, pullback
end

function rrule(::typeof(*), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
    y = Radj * f
    pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
    return y, pullback
end

function rrule(::typeof(\), R::BireStatic, f::Field)
    y = R \ f
    pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
    return y, pullback
end

function rrule(::typeof(\), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
    y = Radj \ f
    pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
    return y, pullback
end


import ChainRulesCore: AbstractTangent

# R * (any tangent)
*(R::BireStatic, t::AbstractTangent) = ZeroTangent()

# R' * (any tangent)
*(Radj::Adjoint{<:Any,<:BireStatic}, t::AbstractTangent) = ZeroTangent()

# R \ (any tangent)
\(R::BireStatic, t::AbstractTangent) = ZeroTangent()

# R' \ (any tangent)
\(Radj::Adjoint{<:Any,<:BireStatic}, t::AbstractTangent) = ZeroTangent()

import ChainRulesCore: rrule, NoTangent, ZeroTangent

function rrule(::Type{<:BireStatic}, α::Field)
    y = BireStatic(α)
    pullback(Δ) = (NoTangent(), ZeroTangent())
    return y, pullback
end

function rrule(::Type{<:BireStatic{T}}, α::Field) where T
    y = BireStatic{T}(α)
    pullback(Δ) = (NoTangent(), ZeroTangent())
    return y, pullback
end