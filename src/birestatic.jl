# abstract type BireStaticOp{T} <: StaticOpWithAdjoint{T} end

# struct BireStatic{T} <: BireStaticOp{T}
#     α::Field
# end

# BireStatic(α::Field) = BireStatic{real(eltype(α))}(α)

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

# # forward, adjoint, inverse and inverse-adjoint
# import Base: *, adjoint, /

# *(R::BireStatic, f::Field) = R(f)
# *(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) = BireStatic(-Radj.parent.α)(f)
# \(R::BireStatic, f::Field) = BireStatic(-R.α)(f)
# \(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) = BireStatic(Radj.parent.α)(f)

# import ChainRulesCore: rrule, NoTangent, ZeroTangent, AbstractTangent

# # function rrule(::typeof(*), R::BireStatic, f::Field)
# #     y = R * f
# #     pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
# #     return y, pullback
# # end

# # function rrule(::typeof(*), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
# #     y = Radj * f
# #     pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
# #     return y, pullback
# # end

# # function rrule(::typeof(\), R::BireStatic, f::Field)
# #     y = R \ f
# #     pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
# #     return y, pullback
# # end

# # function rrule(::typeof(\), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
# #     y = Radj \ f
# #     pullback(Δ) = (NoTangent(), ZeroTangent(), ZeroTangent())
# #     return y, pullback
# # end

# function rrule(::typeof(*), R::BireStatic, f::Field)
#     α = Map(R.α)
#     f_iqu = IQUMap(f)
#     Q = f_iqu.Q
#     U = f_iqu.U
    
#     cos_2α = cos.(2 .* α)
#     sin_2α = sin.(2 .* α)
    
#     Q_rot = cos_2α .* Q .- sin_2α .* U
#     U_rot = sin_2α .* Q .+ cos_2α .* U
    
#     f_iqu_rot = copy(f_iqu)
#     f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
#     f_iqu_rot.arr[:, :, end] .= U_rot.arr
#     y = IQUFourier(f_iqu_rot)
    
#     function pullback(Δ)
#         Δ_iqu = IQUMap(Δ)
#         ΔQ = Δ_iqu.Q
#         ΔU = Δ_iqu.U
        
#         Q_back = cos_2α .* ΔQ .+ sin_2α .* ΔU
#         U_back = -sin_2α .* ΔQ .+ cos_2α .* ΔU
        
#         f_iqu_back = copy(f_iqu)
#         f_iqu_back.arr[:, :, end-1] .= Q_back.arr
#         f_iqu_back.arr[:, :, end] .= U_back.arr
#         ∂f = IQUFourier(f_iqu_back)
        
#         ∂α_from_Q = -2 .* sin_2α .* Q .- 2 .* cos_2α .* U
#         ∂α_from_U = 2 .* cos_2α .* Q .- 2 .* sin_2α .* U
        
#         ∂α = similar(α)
#         ∂α.arr .= ∂α_from_Q.arr .* ΔQ.arr .+ ∂α_from_U.arr .* ΔU.arr
        
#         return (NoTangent(), (α = ∂α,), ∂f)
#     end
    
#     return y, pullback
# end

# function rrule(::typeof(*), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
#     R = Radj.parent
#     α = Map(R.α)
#     f_iqu = IQUMap(f)
#     Q = f_iqu.Q
#     U = f_iqu.U
    
#     cos_2α = cos.(2 .* α)
#     sin_2α = sin.(2 .* α)
    
#     Q_rot = cos_2α .* Q .+ sin_2α .* U
#     U_rot = -sin_2α .* Q .+ cos_2α .* U
    
#     f_iqu_rot = copy(f_iqu)
#     f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
#     f_iqu_rot.arr[:, :, end] .= U_rot.arr
#     y = IQUFourier(f_iqu_rot)
    
#     function pullback(Δ)
#         Δ_iqu = IQUMap(Δ)
#         ΔQ = Δ_iqu.Q
#         ΔU = Δ_iqu.U
        
#         Q_back = cos_2α .* ΔQ .- sin_2α .* ΔU
#         U_back = sin_2α .* ΔQ .+ cos_2α .* ΔU
        
#         f_iqu_back = copy(f_iqu)
#         f_iqu_back.arr[:, :, end-1] .= Q_back.arr
#         f_iqu_back.arr[:, :, end] .= U_back.arr
#         ∂f = IQUFourier(f_iqu_back)
        
#         ∂α_from_Q = 2 .* sin_2α .* Q .- 2 .* cos_2α .* U
#         ∂α_from_U = 2 .* cos_2α .* Q .+ 2 .* sin_2α .* U
        
#         ∂α = similar(α)
#         ∂α.arr .= ∂α_from_Q.arr .* ΔQ.arr .+ ∂α_from_U.arr .* ΔU.arr
        
#         return (NoTangent(), (parent = (α = ∂α,),), ∂f)
#     end
    
#     return y, pullback
# end

# function rrule(::typeof(\), R::BireStatic, f::Field)
#     α = Map(R.α)
#     f_iqu = IQUMap(f)
#     Q = f_iqu.Q
#     U = f_iqu.U
    
#     cos_2α = cos.(2 .* α)
#     sin_2α = sin.(2 .* α)
    
#     Q_rot = cos_2α .* Q .+ sin_2α .* U
#     U_rot = -sin_2α .* Q .+ cos_2α .* U
    
#     f_iqu_rot = copy(f_iqu)
#     f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
#     f_iqu_rot.arr[:, :, end] .= U_rot.arr
#     y = IQUFourier(f_iqu_rot)
    
#     function pullback(Δ)
#         Δ_iqu = IQUMap(Δ)
#         ΔQ = Δ_iqu.Q
#         ΔU = Δ_iqu.U
        
#         Q_back = cos_2α .* ΔQ .- sin_2α .* ΔU
#         U_back = sin_2α .* ΔQ .+ cos_2α .* ΔU
        
#         f_iqu_back = copy(f_iqu)
#         f_iqu_back.arr[:, :, end-1] .= Q_back.arr
#         f_iqu_back.arr[:, :, end] .= U_back.arr
#         ∂f = IQUFourier(f_iqu_back)
        
#         ∂α_from_Q = -2 .* sin_2α .* Q .+ 2 .* cos_2α .* U
#         ∂α_from_U = -2 .* cos_2α .* Q .- 2 .* sin_2α .* U
        
#         ∂α = similar(α)
#         ∂α.arr .= ∂α_from_Q.arr .* ΔQ.arr .+ ∂α_from_U.arr .* ΔU.arr
        
#         return (NoTangent(), (α = ∂α,), ∂f)
#     end
    
#     return y, pullback
# end

# function rrule(::typeof(\), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
#     R = Radj.parent
#     α = Map(R.α)
#     f_iqu = IQUMap(f)
#     Q = f_iqu.Q
#     U = f_iqu.U
    
#     cos_2α = cos.(2 .* α)
#     sin_2α = sin.(2 .* α)
    
#     Q_rot = cos_2α .* Q .- sin_2α .* U
#     U_rot = sin_2α .* Q .+ cos_2α .* U
    
#     f_iqu_rot = copy(f_iqu)
#     f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
#     f_iqu_rot.arr[:, :, end] .= U_rot.arr
#     y = IQUFourier(f_iqu_rot)
    
#     function pullback(Δ)
#         Δ_iqu = IQUMap(Δ)
#         ΔQ = Δ_iqu.Q
#         ΔU = Δ_iqu.U
        
#         Q_back = cos_2α .* ΔQ .+ sin_2α .* ΔU
#         U_back = -sin_2α .* ΔQ .+ cos_2α .* ΔU
        
#         f_iqu_back = copy(f_iqu)
#         f_iqu_back.arr[:, :, end-1] .= Q_back.arr
#         f_iqu_back.arr[:, :, end] .= U_back.arr
#         ∂f = IQUFourier(f_iqu_back)
        
#         ∂α_from_Q = -2 .* sin_2α .* Q .- 2 .* cos_2α .* U
#         ∂α_from_U = 2 .* cos_2α .* Q .- 2 .* sin_2α .* U
        
#         ∂α = similar(α)
#         ∂α.arr .= ∂α_from_Q.arr .* ΔQ.arr .+ ∂α_from_U.arr .* ΔU.arr
        
#         return (NoTangent(), (parent = (α = ∂α,),), ∂f)
#     end
    
#     return y, pullback
# end

# # function rrule(::Type{<:BireStatic}, α::Field)
# #     y = BireStatic(α)
# #     pullback(Δ) = (NoTangent(), ZeroTangent())
# #     return y, pullback
# # end

# # function rrule(::Type{<:BireStatic{T}}, α::Field) where T
# #     y = BireStatic{T}(α)
# #     pullback(Δ) = (NoTangent(), ZeroTangent())
# #     return y, pullback
# # end

# # function rrule(::Type{<:BireStatic}, α::Field)
# #     y = BireStatic(α)
# #     pullback(Δ) = (NoTangent(), Δ.α)
# #     return y, pullback
# # end

# # function rrule(::Type{<:BireStatic{T}}, α::Field) where T
# #     y = BireStatic{T}(α)
# #     pullback(Δ) = (NoTangent(), Δ.α)
# #     return y, pullback
# # end

# using ChainRulesCore: unthunk

# function rrule(::Type{<:BireStatic}, α::Field)
#     y = BireStatic(α)
#     function pullback(Δ)
#         Δ_unwrapped = unthunk(Δ)
#         ∂α = hasproperty(Δ_unwrapped, :α) ? Δ_unwrapped.α : ZeroTangent()
#         return (NoTangent(), ∂α)
#     end
#     return y, pullback
# end

# function rrule(::Type{<:BireStatic{T}}, α::Field) where T
#     y = BireStatic{T}(α)
#     function pullback(Δ)
#         Δ_unwrapped = unthunk(Δ)
#         ∂α = hasproperty(Δ_unwrapped, :α) ? Δ_unwrapped.α : ZeroTangent()
#         return (NoTangent(), ∂α)
#     end
#     return y, pullback
# end

# *(R::BireStatic, t::AbstractTangent) = ZeroTangent()
# *(Radj::Adjoint{<:Any,<:BireStatic}, t::AbstractTangent) = ZeroTangent()
# \(R::BireStatic, t::AbstractTangent) = ZeroTangent()
# \(Radj::Adjoint{<:Any,<:BireStatic}, t::AbstractTangent) = ZeroTangent()

# import ChainRulesCore: AbstractTangent

# *(t::AbstractTangent, ::Any) = ZeroTangent()
# *(::Any, t::AbstractTangent) = ZeroTangent()

# for T in [BireStatic, Adjoint{<:Any,<:BireStatic}]
#     @eval begin
#         *(t::AbstractTangent, ::$T) = ZeroTangent()
#         *( ::$T, t::AbstractTangent) = ZeroTangent()
#     end
# end

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

# forward, adjoint, inverse and inverse-adjoint
import Base: *, adjoint, /

*(R::BireStatic, f::Field) = R(f)
*(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) = BireStatic(-Radj.parent.α)(f)
\(R::BireStatic, f::Field) = BireStatic(-R.α)(f)
\(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) = BireStatic(Radj.parent.α)(f)

import ChainRulesCore: rrule, NoTangent, ZeroTangent, AbstractTangent, Tangent, unthunk

# --- ChainRulesCore rrule definitions ---

function rrule(::typeof(*), R::BireStatic, f::Field)
    α = Map(R.α)
    f_iqu = IQUMap(f)
    Q = f_iqu.Q
    U = f_iqu.U
    
    cos_2α = cos.(2 .* α)
    sin_2α = sin.(2 .* α)
    
    Q_rot = cos_2α .* Q .- sin_2α .* U
    U_rot = sin_2α .* Q .+ cos_2α .* U
    
    f_iqu_rot = copy(f_iqu)
    f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
    f_iqu_rot.arr[:, :, end] .= U_rot.arr
    y = IQUFourier(f_iqu_rot)
    
    function pullback(Δ)
        Δ_iqu = IQUMap(Δ)
        ΔQ = Δ_iqu.Q
        ΔU = Δ_iqu.U
        
        Q_back = cos_2α .* ΔQ .+ sin_2α .* ΔU
        U_back = -sin_2α .* ΔQ .+ cos_2α .* ΔU
        
        f_iqu_back = copy(f_iqu)
        f_iqu_back.arr[:, :, end-1] .= Q_back.arr
        f_iqu_back.arr[:, :, end] .= U_back.arr
        ∂f = IQUFourier(f_iqu_back)
        
        ∂α_from_Q = -2 .* sin_2α .* Q .- 2 .* cos_2α .* U
        ∂α_from_U = 2 .* cos_2α .* Q .- 2 .* sin_2α .* U
        
        ∂α = similar(α)
        ∂α.arr .= ∂α_from_Q.arr .* ΔQ.arr .+ ∂α_from_U.arr .* ΔU.arr
        
        ∂R = Tangent{typeof(R)}(α = ∂α)
        return (NoTangent(), ∂R, ∂f)
    end
    
    return y, pullback
end

function rrule(::typeof(*), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
    R = Radj.parent
    α = Map(R.α)
    f_iqu = IQUMap(f)
    Q = f_iqu.Q
    U = f_iqu.U
    
    cos_2α = cos.(2 .* α)
    sin_2α = sin.(2 .* α)
    
    Q_rot = cos_2α .* Q .+ sin_2α .* U
    U_rot = -sin_2α .* Q .+ cos_2α .* U
    
    f_iqu_rot = copy(f_iqu)
    f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
    f_iqu_rot.arr[:, :, end] .= U_rot.arr
    y = IQUFourier(f_iqu_rot)
    
    function pullback(Δ)
        Δ_iqu = IQUMap(Δ)
        ΔQ = Δ_iqu.Q
        ΔU = Δ_iqu.U
        
        Q_back = cos_2α .* ΔQ .- sin_2α .* ΔU
        U_back = sin_2α .* ΔQ .+ cos_2α .* ΔU
        
        f_iqu_back = copy(f_iqu)
        f_iqu_back.arr[:, :, end-1] .= Q_back.arr
        f_iqu_back.arr[:, :, end] .= U_back.arr
        ∂f = IQUFourier(f_iqu_back)
        
        ∂α_from_Q = 2 .* sin_2α .* Q .- 2 .* cos_2α .* U
        ∂α_from_U = 2 .* cos_2α .* Q .+ 2 .* sin_2α .* U
        
        ∂α = similar(α)
        ∂α.arr .= ∂α_from_Q.arr .* ΔQ.arr .+ ∂α_from_U.arr .* ΔU.arr
        
        ∂R = Tangent{typeof(R)}(α = ∂α)
        return (NoTangent(), ∂R', ∂f)
    end
    
    return y, pullback
end

function rrule(::Type{<:BireStatic}, α::Field)
    y = BireStatic(α)
    function pullback(Δ)
        Δ_unwrapped = unthunk(Δ)
        ∂α = hasproperty(Δ_unwrapped, :α) ? Δ_unwrapped.α : ZeroTangent()
        return (NoTangent(), ∂α)
    end
    return y, pullback
end

function rrule(::Type{<:BireStatic{T}}, α::Field) where T
    y = BireStatic{T}(α)
    function pullback(Δ)
        Δ_unwrapped = unthunk(Δ)
        ∂α = hasproperty(Δ_unwrapped, :α) ? Δ_unwrapped.α : ZeroTangent()
        return (NoTangent(), ∂α)
    end
    return y, pullback
end

*(R::BireStatic, t::AbstractTangent) = ZeroTangent()
*(Radj::Adjoint{<:Any,<:BireStatic}, t::AbstractTangent) = ZeroTangent()
*(t::AbstractTangent, R::BireStatic) = ZeroTangent()
*(t::AbstractTangent, Radj::Adjoint{<:Any,<:BireStatic}) = ZeroTangent()

import ChainRulesCore: AbstractTangent, ZeroTangent

function Base.:*(t::AbstractTangent, Radj::Adjoint{<:Any, <:CachedLenseFlow})
    return ZeroTangent()
end