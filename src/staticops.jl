abstract type StaticOp{T}        <: ImplicitOp{T} end
abstract type StaticOpWithAdjoint{T} <: StaticOp{T} end

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