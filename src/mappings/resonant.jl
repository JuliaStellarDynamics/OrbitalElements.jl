
########################################################################
#
# resonance variables ↔ frequency ratios : Boundaries
#
########################################################################
include("resonant/boundaries.jl")
########################################################################
#
# resonance variables ↔ frequency ratios : mappings
#
########################################################################
"""
    αβ_from_uv(u, v, resonance)

mapping from frequency ratios `(α,β)` to resonance-specific coordinates `(u,v)` given a
resonance number `(n1,n2)`.

See equation (B5) in Fouvry & Prunet (2022)

This signature specifies ωmin and ωmax, to avoid extra calls.
"""
function αβ_from_uv(
    u::Float64,
    v::Float64,            
    res::Resonance
)::Tuple{Float64,Float64}
    n1, n2 = res.number
    ωmin, ωmax = res.frequency_extrema
    if n2 == 0
        return ((ωmax - ωmin) * u + ωmin + ωmax) / (2 * n1), v
    end
    # For v = 0. , β is 0./0.
    # Indeed it implies α = 0 (= αmin) 
    # and u = ±1, ωmin or ωmax = 0., vmin=vmax=0.)
    # In this case, the mapping (α,β) → (u,v) is ill-defined.
    # Every point on the resonance line goes to (ures,0.).
    # We arbritrarly choose to send back all those points on 
    # the radial axis β = 1/2.
    return v, v == 0 ? 1/2 : (((ωmax - ωmin) * u + ωmin + ωmax) / 2 - n1 * v) / (n2 * v)
end

"""
    uv_from_αβ(α, β, resonance)

mapping from resonance-specific coordinates `(u,v)` to frequency ratios `(α,β)` given a
resonance number `(n1,n2)`.

See Appendix B in Fouvry & Prunet (2022)

@IMPROVE, this has rounding error: concern?
"""
function uv_from_αβ(
    α::Float64,
    β::Float64,    
    res::Resonance
)::Tuple{Float64,Float64}
    n1, n2 = res.number
    ωmin, ωmax = res.frequency_extrema
    # Equation B1
    ωval = n1 * α + n2 * α * β
    # Equation B3
    u = (2ωval - ωmax - ωmin) / (ωmax - ωmin)
    # Equation B4
    v = n2 == 0 ? β : α
    return u, v
end

########################################################################
#
# resonance variables ↦ frequency ratios : mapping jacobian
#
########################################################################
"""
    uv_to_αβ_jacobian(u, v, resonance)

jacobian of the (u,v) ↦ (α,β) mapping, i.e. |∂(α,β)/∂(u,v)|
"""
function uv_to_αβ_jacobian(u::Float64,v::Float64,res::Resonance)::Float64
    n1, n2 = res.number
    ωmin, ωmax = res.frequency_extrema
    if n2 == 0
        return (ωmax - ωmin) / abs(2 * n1)
    end        

    return (ωmax - ωmin) / abs(2 * n2 * v)
end

########################################################################
#
# resonance variables ↔ frequency ratios : mappings
#
########################################################################
"""
    rescaled_ϖ(ω̃, resonance)
    rescaled_ϖ(ω̃, ωmin, ωmax)
    rescaled_ϖ(ω̃, n1, n2, model, params)

translate a complex frequency `ω̃` into a resonant rescaled frequency.

See equation (B3) in Fouvry & Prunet (2022)

@ASSUMPTION:
    - `ω̃` is dimensionless, that is, already rescaled by [`frequency_scale(model)`](@ref).
"""
function rescaled_ϖ(
    ω̃::Number,
    n1::Int64,
    n2::Int64,
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)
    ωmin, ωmax = frequency_extrema(n1, n2, model, params)
    return rescaled_ϖ(ω̃, ωmin, ωmax)
end
function rescaled_ϖ(ω̃::Number, res::Resonance)
    ωmin, ωmax = res.frequency_extrema
    return rescaled_ϖ(ω̃, ωmin, ωmax)
end
function rescaled_ϖ(ω̃::Number, ωmin::Float64, ωmax::Float64)
    return (2ω̃ - ωmax - ωmin)/(ωmax - ωmin)
end

########################################################################
#
# resonance lines in action space
#
########################################################################
include("resonant/lines.jl")