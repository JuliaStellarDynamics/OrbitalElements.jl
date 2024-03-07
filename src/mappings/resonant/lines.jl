########################################################################
#
# (J,L) ↦ (u,v) resonances lines in actions space
# These methods for resonance are not use/usable
# @IMPROVE: adapt to v2.0
# @IMPROVE: make them accessible while allowing for other integration methods
# @IMPROVE: Blue style
########################################################################

"""
    actions_resonance_line!(ω, k1, k2, model, tabresult, params)

(k1,k2)-resonance line of the mode in action space.
"""
function actions_resonance_line!(
    ω::Float64,
    k1::Int64,
    k2::Int64,
    model::Potential,
    tabresult::AbstractMatrix{Float64},
    params::OrbitalParameters=OrbitalParameters()
)
    Ω0 = frequency_scale(model)
    # Resonance u
    ωmin, ωmax = Findωminωmax(k1,k2,model,params)
    ures = Getϖ(ω/Ω0,ωmin,ωmax)

    # If no possible resonance (k,k')
    ( -1. < ures < 1.) || (return nothing)

    # Extremal values of v 
    vmin, vmax = FindVminVmax(ures,k1,k2,model,ωmin,ωmax,params)

    # Points along the resonance line
    n, Kv = size(tabresult)
    @assert n == 2 "Wrong result array dimensions"
    δv2 = 1.0/Kv

    for kvval in 1:Kv

        # get the current v2 value
        v2   = δv2 * (kvval - 0.5)
        vval = (vmax - vmin) * (v2^2) + vmin

        ####
        # (u,v) → (J,L)
        ####
        # (u,v) → (α,β)
        α, β = αβ_from_uv(ures, vval, k1, k2, ωmin, ωmax)
        # (α,β) → (Ω1,Ω2)
        Ω1, Ω2 = frequencies_from_αβ(α,β,Ω0)
        # (Ω1,Ω2) → (a,e)
        a, e = ae_from_frequencies(Ω1, Ω2, model, params)
        # (a,e) → (J,L)
        J, L = actions_from_ae(a, e, model, params)

        tabresult[1,kvval], tabresult[2,kvval] = J, L
    end
end

"""
    actions_resonance_lines(ω,model,Kv,tabresonances)

"""
function actions_resonance_lines(
    ω::Number,
    model::CentralPotential,
    Kv::Int64,
    tabresonances::Matrix{Int64},
    params::OrbitalParameters=OrbitalParameters()
)
    # Considered resonances
    n, Kres = size(tabresonances)
    @assert n == 2 "Wrong resonances array dimensions"
    
    # Corresponding actions/frequencies
    tabJLreslines = zeros(Float64,2,Kv,3,Kres)

    ω0, γ = real(ω), imag(ω)

    for ires = 1:Kres

        # Considered resonance pair
        k1, k2   = tabresonances[1,ires], tabresonances[2,ires]

        actions_resonance_line!(ω0-γ,k1,k2,model,view(tabJLreslines,:,:,1,ires),params)
        actions_resonance_line!(ω0,k1,k2,model,view(tabJLreslines,:,:,2,ires),params)
        actions_resonance_line!(ω0+γ,k1,k2,model,view(tabJLreslines,:,:,3,ires),params)
    end

    return tabJLreslines
end