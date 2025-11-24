struct BireFlow
    α::Any
end

function (G::BireFlow)(f)

    E_before = f.E
    B_before = f.B

    α = G.α

    E_new = cos(2 * α) .* E_before .- sin(2 * α) .* B_before
    B_new = sin(2 * α) .* E_before .+ cos(2 * α) .* B_before

    # return E_before, B_before, E_after, B_after

    f_new = copy(f)      # same type, metadata, storage, projection

    # f.arr has shape (Ny/2+1, Nx, 2 [, batch])
    f_new.arr[:, :, 1] .= E_new.arr
    f_new.arr[:, :, 2] .= B_new.arr

    return f_new
end

import Base: *
*(G::BireFlow, f) = G(f)