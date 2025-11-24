abstract type BireFlowOp{T} <: FlowOpWithAdjoint{T} end

struct BireFlow{T} <: BireFlowOp{T}
    α::Field
end

BireFlow(α::Field) = BireFlow{real(eltype(α))}(α)

function (G::BireFlow)(f::Field)
    α = G.α

    f_iqu = IQUFourier(f)

    Q = f_iqu.Q
    U = f_iqu.U

    Q_new = cos.(2 .* α) .* Q .- sin.(2 .* α) .* U
    U_new = sin.(2 .* α) .* Q .+ cos.(2 .* α) .* U

    f_iqu_new = copy(f_iqu)
    f_iqu_new.arr[:, :, end-1] .= Q_new.arr
    f_iqu_new.arr[:, :, end]   .= U_new.arr

    return IEBFourier(f_iqu_new)
end

import Base: *
*(G::BireFlow, f::Field) = G(f)