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

# # Make it at least work, but alpha won't update
# import ChainRulesCore: rrule, NoTangent, ZeroTangent, AbstractTangent, Tangent, unthunk

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
        
#         ∂R = Tangent{typeof(R)}(α = ∂α)
#         return (NoTangent(), ∂R, ∂f)
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
        
#         ∂R = Tangent{typeof(R)}(α = ∂α)
#         return (NoTangent(), ∂R', ∂f)
#     end
    
#     return y, pullback
# end

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

# # Make alpha update correctly
# import Zygote

# Zygote.@adjoint function *(R::BireStatic, f::Field)
#     y = R(f)
    
#     α_field = R.α
#     f_iqu = IQUMap(f)
#     Q = f_iqu.Q
#     U = f_iqu.U
    
#     cos2α = cos.(2 .* α_field)
#     sin2α = sin.(2 .* α_field)
    
#     function back(Δ::Field)
        
#         Δ_iqu = IQUMap(Δ)
#         ΔQ = Δ_iqu.Q 
#         ΔU = Δ_iqu.U
        
#         Q_back = cos2α .* ΔQ .+ sin2α .* ΔU
#         U_back = -sin2α .* ΔQ .+ cos2α .* ΔU
        
#         f_iqu_back = copy(f_iqu)
#         f_iqu_back.arr[:, :, end-1] .= Q_back.arr
#         f_iqu_back.arr[:, :, end]   .= U_back.arr
#         ∇f = IQUFourier(f_iqu_back) 
        
#         M1 = .-sin2α .* ΔQ .+ cos2α .* ΔU
#         M2 = .-cos2α .* ΔQ .- sin2α .* ΔU
        
#         ∇α_field = 2 .* ( Q .* M1 .+ U .* M2 ) 
        
#         ∇R = (α = ∇α_field,)
        
#         return (∇R, ∇f)
#     end
    
#     return y, back
# end

# Zygote.@adjoint function *(Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
#     R = Radj.parent
#     α_field = R.α
    
#     y = BireStatic(-α_field)(f)
    
#     f_iqu = IQUMap(f)
#     Q = f_iqu.Q
#     U = f_iqu.U
    
#     cos2α = cos.(2 .* α_field)
#     sin2α = sin.(2 .* α_field)

#     function back(Δ::Field)
        
#         ∇f = R(Δ)
        
#         Δ_iqu = IQUMap(Δ)
#         ΔQ = Δ_iqu.Q 
#         ΔU = Δ_iqu.U
        
#         M1_neg = sin2α .* ΔQ .+ cos2α .* ΔU
#         M2_neg = .-cos2α .* ΔQ .+ sin2α .* ΔU
        
#         Δneg_α = 2 .* ( Q .* M1_neg .+ U .* M2_neg )
#         ∇α_field = .-Δneg_α
        
#         ∇Radj = (parent = (α = ∇α_field,),)
        
#         return (∇Radj, ∇f)
#     end
    
#     return y, back
# end

# # To surpress certain errors
# import ChainRulesCore: NoTangent, ZeroTangent, AbstractTangent, Tangent

# *(R::BireStatic, t::AbstractTangent) = ZeroTangent()
# *(Radj::Adjoint{<:Any,<:BireStatic}, t::AbstractTangent) = ZeroTangent()
# *(t::AbstractTangent, R::BireStatic) = ZeroTangent()
# *(t::AbstractTangent, Radj::Adjoint{<:Any,<:BireStatic}) = ZeroTangent()

# function Base.:*(t::AbstractTangent, Radj::Adjoint{<:Any, <:CachedLenseFlow})
#     return ZeroTangent()
# end

# function Base.:*(Radj::Adjoint{<:Any, <:CachedLenseFlow}, t::AbstractTangent)
#     return ZeroTangent()
# end

abstract type BireStaticOp{T} <: StaticOpWithAdjoint{T} end

struct BireStatic{T} <: BireStaticOp{T}
    α::Field
end

BireStatic(α::Field) = BireStatic{real(eltype(α))}(α)

function (R::BireStatic)(f::Field)
    α_map = Map(R.α)
    f_iqu = IQUMap(f)

    Q = f_iqu.Q
    U = f_iqu.U

    cos2α = cos.(2 .* α_map)
    sin2α = sin.(2 .* α_map)

    Q_rot = cos2α .* Q .- sin2α .* U
    U_rot = sin2α .* Q .+ cos2α .* U

    f_iqu_rot = copy(f_iqu)
    f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
    f_iqu_rot.arr[:, :, end]   .= U_rot.arr

    return IQUFourier(f_iqu_rot)
end

import Base: *, adjoint, /

*(R::BireStatic, f::Field) = R(f)
*(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) = BireStatic(-Radj.parent.α)(f)
\(R::BireStatic, f::Field) = BireStatic(-R.α)(f)
\(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) = BireStatic(Radj.parent.α)(f)

import ChainRulesCore: rrule, NoTangent, ZeroTangent, AbstractTangent, Tangent, unthunk

function rrule(::typeof(*), R::BireStatic, f::Field)
    α_map = Map(R.α)
    f_iqu = IQUMap(f)
    Q = f_iqu.Q
    U = f_iqu.U

    cos2α = cos.(2 .* α_map)
    sin2α = sin.(2 .* α_map)

    Q_rot = cos2α .* Q .- sin2α .* U
    U_rot = sin2α .* Q .+ cos2α .* U

    f_iqu_rot = copy(f_iqu)
    f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
    f_iqu_rot.arr[:, :, end]   .= U_rot.arr
    y = IQUFourier(f_iqu_rot)

    function pullback(Δ)
        Δ_iqu = IQUMap(Δ)
        ΔQ = Δ_iqu.Q
        ΔU = Δ_iqu.U

        Q_back =  cos2α .* ΔQ .+ sin2α .* ΔU
        U_back = -sin2α .* ΔQ .+ cos2α .* ΔU

        f_iqu_back = copy(f_iqu)
        f_iqu_back.arr[:, :, end-1] .= Q_back.arr
        f_iqu_back.arr[:, :, end]   .= U_back.arr
        ∂f = IQUFourier(f_iqu_back)

        ∂α_from_Q = -2 .* sin2α .* Q .- 2 .* cos2α .* U
        ∂α_from_U =  2 .* cos2α .* Q .- 2 .* sin2α .* U

        ∂α_map = similar(α_map)
        ∂α_map.arr .= ∂α_from_Q.arr .* ΔQ.arr .+ ∂α_from_U.arr .* ΔU.arr

        ∂R = Tangent{typeof(R)}(α = ∂α_map)
        return (NoTangent(), ∂R, ∂f)
    end

    return y, pullback
end

function rrule(::typeof(*), Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
    R = Radj.parent
    α_map = Map(R.α)
    f_iqu = IQUMap(f)
    Q = f_iqu.Q
    U = f_iqu.U

    cos2α = cos.(2 .* α_map)
    sin2α = sin.(2 .* α_map)

    Q_rot =  cos2α .* Q .+ sin2α .* U
    U_rot = -sin2α .* Q .+ cos2α .* U

    f_iqu_rot = copy(f_iqu)
    f_iqu_rot.arr[:, :, end-1] .= Q_rot.arr
    f_iqu_rot.arr[:, :, end]   .= U_rot.arr
    y = IQUFourier(f_iqu_rot)

    function pullback(Δ)
        Δ_iqu = IQUMap(Δ)
        ΔQ = Δ_iqu.Q
        ΔU = Δ_iqu.U

        Q_back = cos2α .* ΔQ .- sin2α .* ΔU
        U_back = sin2α .* ΔQ .+ cos2α .* ΔU

        f_iqu_back = copy(f_iqu)
        f_iqu_back.arr[:, :, end-1] .= Q_back.arr
        f_iqu_back.arr[:, :, end]   .= U_back.arr
        ∂f = IQUFourier(f_iqu_back)

        ∂α_from_Q =  2 .* sin2α .* Q .- 2 .* cos2α .* U
        ∂α_from_U =  2 .* cos2α .* Q .+ 2 .* sin2α .* U

        ∂α_map = similar(α_map)
        ∂α_map.arr .= ∂α_from_Q.arr .* ΔQ.arr .+ ∂α_from_U.arr .* ΔU.arr

        ∂R = Tangent{typeof(R)}(α = ∂α_map)
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

function rrule(::Type{<:BireStatic{T}}, α::Field) where {T}
    y = BireStatic{T}(α)
    function pullback(Δ)
        Δ_unwrapped = unthunk(Δ)
        ∂α = hasproperty(Δ_unwrapped, :α) ? Δ_unwrapped.α : ZeroTangent()
        return (NoTangent(), ∂α)
    end
    return y, pullback
end

import Zygote

Zygote.@adjoint function *(R::BireStatic, f::Field)
    y = R(f)

    α_map = Map(R.α)
    f_iqu = IQUMap(f)
    Q = f_iqu.Q
    U = f_iqu.U

    cos2α = cos.(2 .* α_map)
    sin2α = sin.(2 .* α_map)

    function back(Δ::Field)
        Δ_iqu = IQUMap(Δ)
        ΔQ = Δ_iqu.Q
        ΔU = Δ_iqu.U

        Q_back =  cos2α .* ΔQ .+ sin2α .* ΔU
        U_back = -sin2α .* ΔQ .+ cos2α .* ΔU

        f_iqu_back = copy(f_iqu)
        f_iqu_back.arr[:, :, end-1] .= Q_back.arr
        f_iqu_back.arr[:, :, end]   .= U_back.arr
        ∇f = IQUFourier(f_iqu_back)

        M1 = .-sin2α .* ΔQ .+ cos2α .* ΔU
        M2 = .-cos2α .* ΔQ .- sin2α .* ΔU

        ∇α_map = 2 .* (Q .* M1 .+ U .* M2)

        ∇R = (α = ∇α_map,)
        return (∇R, ∇f)
    end

    return y, back
end

Zygote.@adjoint function *(Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
    R = Radj.parent
    α_map = Map(R.α)

    y = BireStatic(-R.α)(f)

    f_iqu = IQUMap(f)
    Q = f_iqu.Q
    U = f_iqu.U

    cos2α = cos.(2 .* α_map)
    sin2α = sin.(2 .* α_map)

    function back(Δ::Field)
        ∇f = R(Δ)

        Δ_iqu = IQUMap(Δ)
        ΔQ = Δ_iqu.Q
        ΔU = Δ_iqu.U

        M1_neg =  sin2α .* ΔQ .+ cos2α .* ΔU
        M2_neg = .-cos2α .* ΔQ .+ sin2α .* ΔU

        Δneg_α_map = 2 .* (Q .* M1_neg .+ U .* M2_neg)
        ∇α_field = .-Δneg_α_map

        ∇Radj = (parent = (α = ∇α_field,),)
        return (∇Radj, ∇f)
    end

    return y, back
end

# To suppress certain errors
import ChainRulesCore: NoTangent, ZeroTangent, AbstractTangent, Tangent

*(R::BireStatic, t::AbstractTangent) = ZeroTangent()
*(Radj::Adjoint{<:Any,<:BireStatic}, t::AbstractTangent) = ZeroTangent()
*(t::AbstractTangent, R::BireStatic) = ZeroTangent()
*(t::AbstractTangent, Radj::Adjoint{<:Any,<:BireStatic}) = ZeroTangent()

function Base.:*(t::AbstractTangent, Radj::Adjoint{<:Any,<:CachedLenseFlow})
    return ZeroTangent()
end

function Base.:*(Radj::Adjoint{<:Any,<:CachedLenseFlow}, t::AbstractTangent)
    return ZeroTangent()
end