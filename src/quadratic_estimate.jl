export quadratic_estimate, birefringence_quadratic_estimate, joint_quadratic_estimate

#
# Public lensing QE wrapper
#

"""

    quadratic_estimate(ds::DataSet, which; wiener_filtered=true)
    quadratic_estimate((ds₁::DataSet, ds₂::DataSet), which; wiener_filtered=true)

Compute the quadratic estimate of `ϕ` given data.

The `ds` or `(ds₁,ds₂)` tuple contain the DataSet object(s) which house the
data and covariances used in the estimate. Note that only the Fourier-diagonal
approximations for the beam, mask, and noise, i.e. `B̂`, `M̂`, and
`Cn̂`, are accounted for. To account full operators (if they are not actually
Fourier-diagonal), you should compute the impact using Monte Carlo.

If a tuple is passed in, the result will come from correlating the data from
`ds₁` with that from `ds₂`.

An optional keyword argument `AL` can be passed in case the QE normalization
was already computed, in which case it won't be recomputed during the
calculation.

Returns a named tuple of `(;ϕqe, AL, Nϕ)` where `ϕqe` is the (possibly Wiener
filtered, depending on `wiener_filtered` option) quadratic estimate, `AL` is the
normalization (which is already applied to ϕqe, it does not need to be applied
again), and `Nϕ` is the analytic N⁰ noise bias (Nϕ==AL if using unlensed
weights, currently only Nϕ==AL is always returned, no matter the weights)
"""
function quadratic_estimate(
    (ds₁,ds₂) :: NTuple{2,DataSet}, 
    which = nothing; 
    wiener_filtered = true, 
    AL = nothing, 
    weights = :unlensed
)
    @assert weights in (:lensed, :unlensed) "weights should be :lensed or :unlensed"
    if isnothing(which)
        which = ds₁.d isa FlatS0 ? :TT : :EB
    end
    @assert (which in [:TT, :EE, :EB]) "which='$which' not implemented"
    @assert (ds₁.Cf===ds₂.Cf && ds₁.Cf̃===ds₂.Cf̃ && ds₁.Cn̂===ds₂.Cn̂ && ds₁.Cϕ===ds₂.Cϕ && ds₁.B̂===ds₂.B̂) "operators in `ds₁` and `ds₂` should be the same"
    @unpack Cf, Cf̃, Cn̂, Cϕ, B̂, M̂ = ds₁()
    pol = which == :TT ? :I : :P
    quadratic_estimate(Val(which), (ds₁.d[pol], ds₂.d[pol]), Cf[pol], Cf̃[pol], Cn̂[pol], Cϕ, (M̂*B̂)[pol], wiener_filtered, weights, AL)
end

quadratic_estimate(ds::DataSet, args...; kwargs...) = quadratic_estimate((ds,ds), args...; kwargs...)

#
# Birefringence QE
#

"""
    birefringence_quadratic_estimate(ds; weights=:unlensed, zeroB=false)

Compute the unnormalised EB quadratic field for anisotropic birefringence.

For a small polarization rotation,

    Q ± iU -> (Q ± iU) exp(±2iα),

the EB response is

    f_EB^α(l₁,l₂) = 2 [C_EE(l₁) - C_BB(l₂)] cos(2(φ₁-φ₂)).

This function returns the raw, unnormalised quadratic field. Its response,
noise, and possible lensing cross-response should be calibrated by Monte Carlo
or by a separately supplied response matrix.
"""
function birefringence_quadratic_estimate(
    (ds₁, ds₂)::NTuple{2, DataSet};
    weights=:unlensed,
    zeroB=false,
)
    @assert weights in (:lensed, :unlensed) "weights must be :lensed or :unlensed"
    @assert (ds₁.Cf === ds₂.Cf && ds₁.Cf̃ === ds₂.Cf̃ && ds₁.Cn̂ === ds₂.Cn̂ &&
             ds₁.B̂ === ds₂.B̂) "operators in ds₁ and ds₂ must match"

    @unpack Cf, Cf̃, Cn̂, B̂, M̂ = ds₁()
    TF = (M̂ * B̂)[:P]

    birefringence_quadratic_estimate(
        (ds₁.d[:P], ds₂.d[:P]),
        Cf[:P], Cf̃[:P], Cn̂[:P], TF;
        weights,
        zeroB,
    )
end

birefringence_quadratic_estimate(ds::DataSet; kwargs...) =
    birefringence_quadratic_estimate((ds, ds); kwargs...)


function birefringence_quadratic_estimate(
    (d₁, d₂)::NTuple{2, FlatS2},
    Cf, Cf̃, Cn, TF;
    weights=:unlensed,
    zeroB=false,
)
    @assert weights in (:lensed, :unlensed) "weights must be :lensed or :unlensed"

    Cw = (weights == :unlensed) ? Cf : Cf̃
    CE, CB = Cw[:E], Cw[:B]

    TF²E, TF²B = TF[:E]^2, TF[:B]^2
    ΣEtot = TF²E * Cf̃[:E] + Cn[:E]
    ΣBtot = TF²B * Cf̃[:B] + Cn[:B]

    E₁ = (ΣEtot \ (TF * d₁)[:E])
    B₂ = (ΣBtot \ (TF * d₂)[:B])

    # First response term:
    #   + 2 C_EE(l_E) cos(2Δφ) E(l_E) B(l_B)
    E_CE = CE * E₁
    B_id = B₂

    cE_CE = QE_leg(E_CE, 1, 1) - QE_leg(E_CE, 2, 2)
    sE_CE = 2 * QE_leg(E_CE, 1, 2)

    cB_id = QE_leg(B_id, 1, 1) - QE_leg(B_id, 2, 2)
    sB_id = 2 * QE_leg(B_id, 1, 2)

    term_EE = cE_CE * cB_id + sE_CE * sB_id

    # Second response term:
    #   - 2 C_BB(l_B) cos(2Δφ) E(l_E) B(l_B)
    if zeroB
        term_BB = 0
    else
        E_id = E₁
        B_CB = CB * B₂

        cE_id = QE_leg(E_id, 1, 1) - QE_leg(E_id, 2, 2)
        sE_id = 2 * QE_leg(E_id, 1, 2)

        cB_CB = QE_leg(B_CB, 1, 1) - QE_leg(B_CB, 2, 2)
        sB_CB = 2 * QE_leg(B_CB, 1, 2)

        term_BB = cE_id * cB_CB + sE_id * sB_CB
    end

    # Full EB birefringence raw field:
    #
    #   α_raw ∝ 2 [C_EE(l_E) - C_BB(l_B)] cos(2Δφ) E B
    #
    αqe = 2 * Fourier(term_EE - term_BB)

    Memoization.empty_cache!(QE_leg)

    (; αqe)
end

#
# QE leg helpers
#

@doc doc"""

    QE_leg(C::Diagonal, inds...)

The quadratic estimate and normalization expressions all consist of terms
involving products of two "legs", each of the form:

    C * l[i] * l̂[j] * l̂[k] * ...

where `C` is a field or diagonal covariance, `l[i]` is the Fourier wave-vector
in direction `i`, and `l̂[i] = l[i] / ‖l‖`. For example, an EB-estimator leg
like

    (CE * (CẼ + Cn) \ d[:E]) * l[i] * l̂[j] * l̂[k]

is represented as:

    QE_leg((CE * (CẼ + Cn) \ d[:E]), [i], j, k)

Putting an index in brackets selects the Fourier wave-vector `l`; a plain
integer selects the unit-vector `l̂`.

The legs are symmetric in their indices. `QE_leg` memoizes each unique index
structure, which avoids repeated work in EE and EB normalization terms.
"""
QE_leg(C::Diagonal, inds...) = QE_leg(C.diag, inds...)
function QE_leg(C::Field, inds...)
    n = count((x -> x isa Int), inds)
    p₁, p₂ = (count(==(i), first.(inds)) for i = 1:2)
    QE_leg(C, (n, p₁, p₂))
end
@memoize function QE_leg(C::Field, (n, p₁, p₂)::Tuple)
    Map(@. nan2zero($Ð(C) * ∇[1].diag^p₁ * ∇[2].diag^p₂ / sqrt(∇².diag)^n))
end
ϵ(x...) = levicivita([x...])
inds(D) = collect(product(repeated(1:2, D)...))[:]

#
# Lensing QE implementations
#

function quadratic_estimate(::Val{:TT}, (d₁, d₂)::NTuple{2, FlatS0}, Cf, Cf̃, Cn, Cϕ, TF, wiener_filtered, weights, AL=nothing)

    ΣTtot = TF^2 * Cf̃ + Cn
    CT = (weights == :unlensed) ? Cf : Cf̃

    ϕqe_unnormalized = @subst -sum(∇[i] * Fourier(QE_leg($(ΣTtot \ (TF * d₁))) * QE_leg($(CT * (ΣTtot \ (TF * d₂))), [i])) for i = 1:2)

    if AL == nothing
        AL = @subst begin
            A(i, j) = (
                QE_leg($(TF^2 * CT^2 / ΣTtot), [i], [j]) * QE_leg($(TF^2 / ΣTtot))
              + QE_leg($(TF^2 * CT / ΣTtot), [i]) * QE_leg($(TF^2 * CT / ΣTtot), [j])
            )
            pinv(Diagonal(sum(abs.(∇[i].diag .* ∇[j].diag .* Fourier(A(i, j))) for (i, j) in inds(2))))
        end
    end
    Nϕ = AL

    Memoization.empty_cache!(QE_leg)

    ϕqe = (wiener_filtered ? (Cϕ * pinv(Cϕ + Nϕ)) : 1) * (AL * ϕqe_unnormalized)
    (; ϕqe, AL, Nϕ)

end

function quadratic_estimate(::Val{:EE}, (d₁, d₂)::NTuple{2, FlatS2}, Cf, Cf̃, Cn, Cϕ, TF, wiener_filtered, weights, AL=nothing)

    TF² = TF[:E]^2
    ΣEtot = TF² * Cf̃[:E] + Cn[:E]
    CE = ((weights == :unlensed) ? Cf : Cf̃)[:E]

    ϕqe_unnormalized = @subst begin
        I(i) = -(
            2sum(QE_leg($(CE * (ΣEtot \ (TF * d₁)[:E])), [i], j, k) * QE_leg($((ΣEtot \ (TF * d₂)[:E])), j, k) for (j, k) in inds(2))
               - QE_leg($(CE * (ΣEtot \ (TF * d₁)[:E])), [i]) * QE_leg($((ΣEtot \ (TF * d₂)[:E])))
        )
        sum(∇[i] * Fourier(I(i)) for i = 1:2)
    end

    if AL == nothing
        AL = @subst begin
            A1(i, j) = -4 * sum(ϵ(m, p, 3) * ϵ(n, q, 3) * (
                  QE_leg($(TF² * CE^2 / ΣEtot), [i], [j], k, l, m, n) .* QE_leg($(TF² / ΣEtot), k, l, p, q)
                + QE_leg($(TF² * CE / ΣEtot), [i], k, l, m, n) .* QE_leg($(TF² * CE / ΣEtot), [j], k, l, p, q))
                for (k, l, m, n, p, q) in inds(6)
            )
            A2(i, j) = (
                  QE_leg($(TF² * CE^2 / ΣEtot), [i], [j]) .* QE_leg($(TF² / ΣEtot))
                + QE_leg($(TF² * CE / ΣEtot), [i]) .* QE_leg($(TF² * CE / ΣEtot), [j])
            )
            pinv(Diagonal(sum(abs.(∇[i].diag .* ∇[j].diag .* Fourier(A1(i, j) + A2(i, j))) for (i, j) in inds(2))))
        end
    end
    Nϕ = AL

    Memoization.empty_cache!(QE_leg)

    ϕqe = (wiener_filtered ? (Cϕ * pinv(Cϕ + Nϕ)) : 1) * (AL * ϕqe_unnormalized)
    (; ϕqe, AL, Nϕ)

end

function quadratic_estimate(::Val{:EB}, (d₁, d₂)::NTuple{2, FlatS2}, Cf, Cf̃, Cn, Cϕ, TF, wiener_filtered, weights, AL=nothing; zeroB=false)

    CE, CB = getindex.(Ref((weights == :unlensed) ? Cf : Cf̃), (:E, :B))
    TF²E, TF²B = TF[:E]^2, TF[:B]^2
    ΣEtot = TF²E * Cf̃[:E] + Cn[:E]
    ΣBtot = TF²B * Cf̃[:B] + Cn[:B]

    ϕqe_unnormalized = @subst begin
        I(i) = 2 * sum(ϵ(k, l, 3) * (
                           QE_leg($(CE * (ΣEtot \ (TF * d₁)[:E])), [i], j, k) * QE_leg($((ΣBtot \ (TF * d₂)[:B])), j, l)
            - (zeroB ? 0 : QE_leg($(ΣEtot \ (TF * d₁)[:E]), j, k) * QE_leg($(CB * (ΣBtot \ (TF * d₂)[:B])), [i], j, l)))
            for (j, k, l) in inds(3)
        )
        sum(∇[i] * Fourier(I(i)) for i = 1:2)
    end

    if AL == nothing
        AL = @subst begin
            A(i, j) = 4 * sum(ϵ(m, p, 3) * ϵ(n, q, 3) * (
                                 QE_leg($(TF²E * CE^2 / ΣEtot), [i], [j], k, l, m, n) * QE_leg($(TF²B / ΣBtot), k, l, p, q)
                + (zeroB ? 0 : -2QE_leg($(TF²E * CE / ΣEtot), [i], k, l, m, n) * QE_leg($(TF²B * CB / ΣBtot), [j], k, l, p, q))
                + (zeroB ? 0 :   QE_leg($(TF²E / ΣEtot), k, l, m, n) * QE_leg($(TF²B * CB^2 / ΣBtot), [i], [j], k, l, p, q)))
                for (k, l, m, n, p, q) in inds(6)
            )
            pinv(Diagonal(sum(abs.(∇[i].diag .* ∇[j].diag .* Fourier(A(i, j))) for (i, j) in inds(2))))
        end
    end
    Nϕ = AL

    Memoization.empty_cache!(QE_leg)

    ϕqe = (wiener_filtered ? (Cϕ * pinv(Cϕ + Nϕ)) : 1) * (AL * ϕqe_unnormalized)
    (; ϕqe, AL, Nϕ)

end

#
# Joint lensing and birefringence wrapper
#

"""
    joint_quadratic_estimate(ds;
        phi_weights=:unlensed,
        alpha_weights=:unlensed,
        wiener_filtered=false,
        zeroB_alpha=false,
        response=nothing,
    )

Return EB lensing and birefringence quadratic fields from the same data.

If `response === nothing`, this returns the separate outputs:

    ϕ = quadratic_estimate(ds, :EB; ...)
    α = birefringence_quadratic_estimate(ds; ...)

If `response` is supplied, it should be a named tuple

    (; Rϕϕ, Rϕα, Rαϕ, Rαα)

representing the MC-calibrated 2x2 response matrix in each L mode:

    [ϕ_raw]   [Rϕϕ  Rϕα] [ϕ]
    [α_raw] = [Rαϕ  Rαα] [α].

The returned `ϕhat` and `αhat` are response-matrix-deprojected estimates.
"""
function joint_quadratic_estimate(
    ds::DataSet;
    phi_weights=:unlensed,
    alpha_weights=:unlensed,
    wiener_filtered=false,
    zeroB_alpha=false,
    response=nothing,
)
    ϕ = quadratic_estimate(
        ds, :EB;
        weights=phi_weights,
        wiener_filtered=wiener_filtered,
    )

    α = birefringence_quadratic_estimate(
        ds;
        weights=alpha_weights,
        zeroB=zeroB_alpha,
    )

    ϕraw = ϕ.ϕqe
    αraw = α.αqe

    if response === nothing
        return (; ϕ, α)
    end

    @assert all(k -> haskey(response, k), (:Rϕϕ, :Rϕα, :Rαϕ, :Rαα)) """
    response must be a named tuple with fields:
        (; Rϕϕ, Rϕα, Rαϕ, Rαα)
    """

    Rϕϕ = response.Rϕϕ
    Rϕα = response.Rϕα
    Rαϕ = response.Rαϕ
    Rαα = response.Rαα

    detR = Rϕϕ * Rαα - Rϕα * Rαϕ

    ϕhat = ( Rαα * ϕraw - Rϕα * αraw) / detR
    αhat = (-Rαϕ * ϕraw + Rϕϕ * αraw) / detR

    (; ϕhat, αhat, ϕraw, αraw, ϕ, α, response)
end
