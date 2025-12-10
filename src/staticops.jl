# abstract type StaticOp{T}        <: ImplicitOp{T} end
# abstract type StaticOpWithAdjoint{T} <: StaticOp{T} end

# *(A::StaticOp, f::Field) = A(f)

# *(Aadj::Adjoint{<:Any,<:StaticOp}, f::Field) =
#     error("Adjoint not implemented for $(typeof(Aadj.parent))")

# \(A::StaticOp, f::Field) =
#     error("Inverse not implemented for $(typeof(A))")

# \(Aadj::Adjoint{<:Any,<:StaticOp}, f::Field) =
#     error("Inverse-adjoint not implemented for $(typeof(Aadj.parent))")

# import Zygote

# Zygote.@adjoint function *(A::StaticOp, f::Field)
#     y = A(f)
#     function back(Δ)
#         # derivative wrt operator is nothing
#         # derivative wrt f is A' * Δ
#         return (nothing, A' * Δ)
#     end
#     return y, back
# end

# Zygote.@adjoint function \(A::StaticOp, f::Field)
#     y = A \ f
#     function back(Δ)
#         # derivative wrt operator is nothing
#         # derivative wrt f is A' \ Δ
#         return (nothing, A' \ Δ)
#     end
#     return y, back
# end

abstract type StaticOp{T}        <: ImplicitOp{T} end
abstract type StaticOpWithAdjoint{T} <: StaticOp{T} end

*(A::StaticOp, f::Field) = A(f)

*(Aadj::Adjoint{<:Any,<:StaticOp}, f::Field) =
    error("Adjoint not implemented for $(typeof(Aadj.parent))")

\(A::StaticOp, f::Field) =
    error("Inverse not implemented for $(typeof(A))")

\(Aadj::Adjoint{<:Any,<:StaticOp}, f::Field) =
    error("Inverse-adjoint not implemented for $(typeof(Aadj.parent))")

import Zygote

Zygote.@adjoint function *(A::StaticOp, f::Field)
    y = A(f)
    function back(Δ)
        return (nothing, A' * Δ)
    end
    return y, back
end

Zygote.@adjoint function \(A::StaticOp, f::Field)
    y = A \ f
    function back(Δ)
        return (nothing, A' \ Δ)
    end
    return y, back
end

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

import Base: *, adjoint, /

*(R::BireStatic, f::Field) = R(f)
*(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) = BireStatic(-Radj.parent.α)(f)
\(R::BireStatic, f::Field) = BireStatic(-R.α)(f)
\(Radj::Adjoint{<:Any,<:BireStatic}, f::Field) = BireStatic(Radj.parent.α)(f)

Zygote.@adjoint function *(R::BireStatic, f::Field)
    y = R(f)
    
    α_field = R.α
    f_iqu = IQUMap(f)
    Q = f_iqu.Q
    U = f_iqu.U
    
    cos2α = cos.(2 .* α_field)
    sin2α = sin.(2 .* α_field)
    
    function back(Δ::Field)
        
        Δ_iqu = IQUMap(Δ)
        ΔQ = Δ_iqu.Q 
        ΔU = Δ_iqu.U
        
        Q_back = cos2α .* ΔQ .+ sin2α .* ΔU
        U_back = -sin2α .* ΔQ .+ cos2α .* ΔU
        
        f_iqu_back = copy(f_iqu)
        f_iqu_back.arr[:, :, end-1] .= Q_back.arr
        f_iqu_back.arr[:, :, end]   .= U_back.arr
        ∇f = IQUFourier(f_iqu_back) 
        
        M1 = .-sin2α .* ΔQ .+ cos2α .* ΔU
        M2 = .-cos2α .* ΔQ .- sin2α .* ΔU
        
        ∇α_field = 2 .* ( Q .* M1 .+ U .* M2 ) 
        
        ∇R = (α = ∇α_field,)
        
        return (∇R, ∇f)
    end
    
    return y, back
end

Zygote.@adjoint function *(Radj::Adjoint{<:Any,<:BireStatic}, f::Field)
    R = Radj.parent
    α_field = R.α
    
    y = BireStatic(-α_field)(f)
    
    f_iqu = IQUMap(f)
    Q = f_iqu.Q
    U = f_iqu.U
    
    cos2α = cos.(2 .* α_field)
    sin2α = sin.(2 .* α_field)

    function back(Δ::Field)
        
        ∇f = R(Δ)
        
        Δ_iqu = IQUMap(Δ)
        ΔQ = Δ_iqu.Q 
        ΔU = Δ_iqu.U
        
        M1_neg = sin2α .* ΔQ .+ cos2α .* ΔU
        M2_neg = .-cos2α .* ΔQ .+ sin2α .* ΔU
        
        Δneg_α = 2 .* ( Q .* M1_neg .+ U .* M2_neg )
        ∇α_field = .-Δneg_α
        
        ∇Radj = (parent = (α = ∇α_field,),)
        
        return (∇Radj, ∇f)
    end
    
    return y, back
end