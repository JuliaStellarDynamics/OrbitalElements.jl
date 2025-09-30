########################################################################
#
# (J,L) ↦ (u,v) resonances lines in actions space
# These methods for resonance are not use/usable
# @WARNING: This file is not included in the main code anyway !
# @IMPROVE: adapt to v2.0
# @IMPROVE: make them accessible while allowing for other integration methods
# @IMPROVE: Blue style
########################################################################

########################################################################
#
# resonance variables : rescaling v into v'
#
########################################################################
const DEFAULT_VMAPN = 2
"""
    v_from_vp(vp, vmin, vmax[, n])

resonant variable `v` as a function of its dilated/rescaled version `vp`
"""
function v_from_vp(vp::Float64, vmin::Float64, vmax::Float64; n::Int64=DEFAULT_VMAPN)
    return (vmax - vmin) * vp^n + vmin
end

"""
    v_from_vp_derivative(vp, vmin, vmax[, n])

same as `v_from_vp` plus mapping derivative
"""
function v_from_vp_derivative(
    vp::Float64, vmin::Float64, vmax::Float64; n::Int64=DEFAULT_VMAPN
) 
    return v_from_vp(vp, vmin, vmax, n=n), n * (vmax - vmin) * vp^(n-1)
end

"""
    v_from_vp(vp, vmin, vmax[, n])

resonant dilated/rescaled `vp` variable as a function of `v`
"""
function vp_from_v(v::Float64, vmin::Float64, vmax::Float64; n::Int64=DEFAULT_VMAPN)
    return ((v - vmin) / (vmax - vmin))^(1/n)
end

"""
    actions_resonance_line!(result, ω, resonance, model, params[, vmapn])

resonance line of  in action space, in place.
"""
function actions_resonance_line!(
    result::AbstractMatrix{Float64},
    ω::Float64,
    resonance::Resonance,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters();
    vmapn::Int64=DEFAULT_VMAPN
)
    # Check result array of dimension (2, Kv)
    dimJ, Kv = size(result)
    @assert dimJ == 2 "Wrong result array dimensions"

    fill!(result, 0.0)
    # ω has to be rescaled by the model frequency
    ω̃ = ω / frequency_scale(model)
    ures = rescaled_ϖ(ω̃, resonance)

    # If no possible resonance
    ( -1 < ures < 1) || (return result)

    # Extremal values of v 
    vmin, vmax = v_boundaries(ures, resonance, model, params)

    # Number of points along the resonance line
    δvp = 1 / Kv

    for kvval in 1:Kv

        # get the current v2 value
        vp = δvp * (kvval - 0.5)
        vval = v_from_vp(vp, vmin, vmax, n=vmapn)

        ####
        # (u,v) → (J,L)
        ####
        # (u,v) → (α,β)
        α, β = αβ_from_uv(ures, vval, resonance)
        # (α,β) → (Ω1,Ω2)
        Ω1, Ω2 = frequencies_from_αβ(α, β, model)
        # (Ω1,Ω2) → (a,e)
        a, e = ae_from_frequencies(Ω1, Ω2, model, params)
        # (a,e) → (J,L)
        J, L = actions_from_ae(a, e, model, params)

        result[1,kvval], result[2,kvval] = J, L
    end

    return result
end

const DEFAULT_KV = 100
"""
    actions_resonance_line(ω, resonance, model, params[, Kv, vmapn])

same as `actions_resonance_line!` but allocating.
"""
function actions_resonance_line(
    ω::Float64,
    resonance::Resonance,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters();
    Kv::Int64=DEFAULT_KV,
    vmapn::Int64=DEFAULT_VMAPN,
)
    # allocate
    result = zeros(2, Kv)
    return actions_resonance_line!(result, ω, resonance, model, params, vmapn=vmapn)
end