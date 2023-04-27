
"""
    GetResLineJL!(ω,k1,k2,ψ,dψ,d2ψ,tabresult,params)

(k1,k2)-resonance line of the mode in action space.
"""
function GetResLineJL!(ω::Float64,
                       k1::Int64,k2::Int64,
                       ψ::F0,dψ::F1,d2ψ::F2,
                       tabresult::AbstractMatrix{Float64},
                       params::OrbitalParameters=OrbitalParameters()) where {F0 <: Function, F1 <: Function, F2 <: Function}

    Ω₀ = params.Ω₀

    # Resonance u
    ωmin, ωmax = Findωminωmax(k1,k2,dψ,d2ψ,params)
    ures = Getϖ(ω/Ω₀,ωmin,ωmax)

    # If no possible resonance (k,k')
    ( -1. < ures < 1.) || (return nothing)

    # Extremal values of v 
    vmin, vmax = FindVminVmax(ures,k1,k2,dψ,d2ψ,ωmin,ωmax,params)

    # Points along the resonance line
    n, Kv = size(tabresult)
    @assert n == 2 "Wrong result array dimensions"
    δv2 = 1.0/Kv

    for kvval in 1:Kv

        # get the current v2 value
        v2   = δv2*(kvval-0.5)
        vval = (vmax-vmin)*(v2^2)+vmin

        ####
        # (u,v) → (J,L)
        ####
        # (u,v) → (α,β)
        α, β = αβFromUV(ures,vval,k1,k2,ωmin,ωmax)
        # (α,β) → (Ω1,Ω2)
        Ω1, Ω2 = FrequenciesFromαβ(α,β,Ω₀)
        # (Ω1,Ω2) → (a,e)
        a, e = ComputeAEFromFrequencies(ψ,dψ,d2ψ,Ω1,Ω2,params)
        # (a,e) → (J,L)
        J, L = ComputeActionsAE(ψ,dψ,a,e,params)

        tabresult[1,kvval], tabresult[2,kvval] = J, L
    end
end

"""
    GetResLinesJL(ω,ψ,dψ,d2ψ,Kv,tabresonances)

"""
function GetResLinesJL(ω::Number,
                       ψ::F0,dψ::F1,d2ψ::F2,
                       Kv::Int64,tabresonances::Matrix{Int64},
                       params::OrbitalParameters=OrbitalParameters()) where {F0 <: Function, F1 <: Function, F2 <: Function}

    # Considered resonances
    n, Kres = size(tabresonances)
    @assert n == 2 "Wrong resonances array dimensions"
    
    # Corresponding actions/frequencies
    tabJLreslines = zeros(Float64,2,Kv,3,Kres)

    ω0, γ = real(ω), imag(ω)

    for ires = 1:Kres

        # Considered resonance pair
        k1, k2   = tabresonances[1,ires], tabresonances[2,ires]

        GetResLineJL!(ω0-γ,k1,k2,ψ,dψ,d2ψ,view(tabJLreslines,:,:,1,ires),params)
        GetResLineJL!(ω0,k1,k2,ψ,dψ,d2ψ,view(tabJLreslines,:,:,2,ires),params)
        GetResLineJL!(ω0+γ,k1,k2,ψ,dψ,d2ψ,view(tabJLreslines,:,:,3,ires),params)
    end

    return tabJLreslines
end