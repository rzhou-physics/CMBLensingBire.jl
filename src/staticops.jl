# ============================================
#  Static linear operators
#
#  These apply a linear map f ↦ A f without 
#  time evolution (unlike FlowOp).
#  Examples: Beam, Mask, Filter, BireStatic
# ============================================

abstract type StaticOp{T}        <: ImplicitOp{T} end
abstract type StaticOpWithAdjoint{T} <: StaticOp{T} end

# --------------------------------------------
#  Forward operator:  A * f
# --------------------------------------------
*(A::StaticOp, f::Field) = A(f)


# --------------------------------------------
#  Adjoint operator:  A' * f
#
#  The user must define:
#     *(Adjoint(A), f)
# --------------------------------------------
*(Aadj::Adjoint{<:Any,<:StaticOp}, f::Field) =
    error("Adjoint not implemented for $(typeof(Aadj.parent))")


# --------------------------------------------
#  Inverse operator:  A \ f
#
#  The user must implement:
#     \(A::StaticOp, f::Field)
# --------------------------------------------
\(A::StaticOp, f::Field) =
    error("Inverse not implemented for $(typeof(A))")


# --------------------------------------------
#  Inverse-adjoint operator:  A' \ f
#
#  The user must implement:
#     \(Adjoint(A), f)
# --------------------------------------------
\(Aadj::Adjoint{<:Any,<:StaticOp}, f::Field) =
    error("Inverse-adjoint not implemented for $(typeof(Aadj.parent))")


# --------------------------------------------
#  Automatic differentiation w.r.t. f
#
#  NOTE: StaticOpWithAdjoint *provides its own adjoint* 
#        for the operator argument.
#
#  We only differentiate w.r.t. f here.
# --------------------------------------------

import Zygote

# d/d f of A * f
Zygote.@adjoint function *(A::StaticOp, f::Field)
    y = A(f)
    function back(Δ)
        # derivative wrt operator is nothing
        # derivative wrt f is A' * Δ
        return (nothing, A' * Δ)
    end
    return y, back
end

# d/d f of A \ f
Zygote.@adjoint function \(A::StaticOp, f::Field)
    y = A \ f
    function back(Δ)
        # derivative wrt operator is nothing
        # derivative wrt f is A' \ Δ
        return (nothing, A' \ Δ)
    end
    return y, back
end