
########################################################################
#
# Appropriate reduced coordinates
#
# @IMPROVE: Merge with Isochrone !! Remove code duplicates
#           + Confusing names among anomalies: use s for this specific anomaly ?      
#
#
########################################################################
"""
    _anomaly_from_radius(r, model::PlummerPotential)

anomaly independent of the orbit (only radius dependent) for Plummer.

@IMPROVE: confusing names ! anomaly `s` vs Henon anomaly `u`.
"""
function _anomaly_from_radius(r::Float64, model::PlummerPotential)
    x = r / radial_scale(model) # dimensionless radius
    return sqrt(1 + x^2)
end

"""
    _radius_from_anomaly(anomaly, model::PlummerPotential)

anomaly independent of the orbit (only radius dependent) for Plummer.

@IMPROVE: confusing names ! anomaly `s` vs Henon anomaly `u`.
"""
function _radius_from_anomaly(anomaly::Float64, model::SemiAnalyticPlummer)
    x = sqrt(anomaly^2 - 1) # dimensionless radius
    return radial_scale(model) * x
end

"""
    _spsa_from_rpra(rp, ra, model::PlummerPotential)

extremal anomalies of the orbit with pericentre `rp` and apocentre `ra`.
"""
function _spsa_from_rpra(rp::Float64, ra::Float64, model::SemiAnalyticPlummer)
    # Broadcast the anomaly over peri and apocenter (allocation-free)
    sp, sa = map((r -> _anomaly_from_radius(r, model)), (rp, ra))
    return sp, sa
end

"""
    _rpra_from_spsa(sp, sa, model::PlummerPotential)

the orbit pericentre `rp` and apocentre `ra` from extremal anomalies.
"""
function _rpra_from_spsa(sp::Float64, sa::Float64, model::SemiAnalyticPlummer)
    # Broadcast the anomaly over peri and apocenter (allocation-free)
    rp, ra = map((s -> _radius_from_anomaly(s, model)), (sp, sa))
    return rp, ra
end

########################################################################
#
# (a,e) ↔ (E,L) mappings
#
########################################################################
"""
for Plummer analytical version, see equation @ADDREFERENCE
"""
function EL_from_ae(
    a::Float64,
    e::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)
    rp, ra = rpra_from_ae(a, e)
    sp, sa = _spsa_from_rpra(rp, ra, model) # extremal anomalies on the orbit
    Ẽ = 1 / sp - (sa^2 - 1) / (sa * sp * (sa + sp)) # dimensionless energy
    L̃ = sqrt(2 * (sp^2 - 1) * (sa^2 - 1) / (sa * sp * (sa + sp))) # dimensionless momentum
    return energy_scale(model) * Ẽ, momentum_scale(model) * L̃ 
end

"""    
for Plummer semi-analytical version, can be found by 1D-bisection instead of 2D inversion scheme.
"""
function ae_from_EL(
    E::Float64,
    L::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)
    sp, sa = _spsa_from_EL(E, L, model, params)
    rp, ra = _rpra_from_spsa(sp, sa, model)
    return ae_from_rpra(rp, ra)
end

function _spsa_from_EL(
    E::Float64,
    L::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)
    Ẽ = E / energy_scale(model) # dimensionless energy
    L̃ = L / momentum_scale(model) # dimensionless angular momentum

    # Edge cases
    # @ADDREFERENCE for this cryptic formula to find root from
    rootequation = (s -> Ẽ * s^3 - s^2 + ((L̃^2) / 2 - Ẽ) * s + 1)
    if Ẽ <= 0 # unbounded orbits
        if L̃ == 0 # radial orbits
            return 1.0, Inf
        else # arbitrary orbit
            sp = _bisection(rootequation, 1.0, 2 / (L̃^2 - 2Ẽ))
            return sp, Inf
        end
    elseif L̃ == 0 # radial orbit
        return 1.0, 1 / Ẽ
    end

    #  Generic case: bounded orbits
    # 1 < sp < sc and sc < sa < E0/E with sc the circular anomaly
    # @ADDREFERENCE for this cryptic formulae
    # Anomaly of the circular orbit of energy E
    x1, x2, x3 = 1 / (6Ẽ), (1 + 54Ẽ^2) / (216Ẽ^3), sqrt(1 + 1 / (27Ẽ^2)) / (4Ẽ)
    scircular = x1 + ∛(x2 + x3) + ∛(x2 - x3)
    η̃ = (scircular^2 - 1) * (1 / scircular - Ẽ)
    L̃circular = Ẽ == 1 ? 0.0 : sqrt(abs(2η̃))
    if L̃ >= L̃circular # circular orbits
        # @IMPROVE: inequality to take care of small Float errors, not ideal
        # should find some proper cutoff for quasi circular orbits
        return scircular, scircular
    end

    sp = _bisection(rootequation, 1.0, scircular)
    # upper bound should be E0/E + 1 in order to get a proper bracket
    sa = _bisection(rootequation, scircular, 1 + 1/Ẽ)
    return sp, sa
end

########################################################################
#
# (a,e) ↔ frequencies mappings: integrands
#
########################################################################
"""
the analytic expression for Theta for the Plummer profile are given in 
Tep et al. 2022 eq. (F9). It has the major advantage of 
always being well-posed for -1<u<1

@QUESTION: is `u` hard-coded as Hénon's anomaly here ?
"""
function Θ(
    u::Float64,
    a::Float64,
    e::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)
    rp, ra = rpra_from_ae(a, e) # peri and apocentre
    sp, sa = _spsa_from_rpra(rp, ra, model) # extremal anomalies
    
    # Analytical expression of (dr/du)(1/vr), that is always well-posed
    # see equation (F9) in Tep et al. (2022).
    tep_a = sp * (u + 2) * (u - 1)^2 - sa * (u - 2) * (u + 1)^2
    tep_b = sp * (u^3 - 3u + 6) - sa * (u^3 - 3u - 6)

    return (
        3 
        * sqrt(sa * sp * (sa + sp))
        * tep_a^(3/2)
        / sqrt(sa * sp * tep_a + tep_b)
        / sqrt(4 - u^2)
        / (2 * sqrt(2) * frequency_scale(model))
    )
end

########################################################################
#
# (a, e) ↔ frequencies mappings: circular
#
########################################################################
function _Ω1circular(a::Float64, model::SemiAnalyticPlummer)::Float64
    x = a / radial_scale(model) # Dimensionless radius
    anomaly = _anomaly_from_radius(a, model)
    return frequency_scale(model) * sqrt(4 + x^2) / (2 * anomaly^(5/2))
end

function _Ω2circular(a::Float64, model::SemiAnalyticPlummer)::Float64
    anomaly = _anomaly_from_radius(a, model)
    return frequency_scale(model) / (2 * anomaly^(3/2))
end

function _βcircular(a::Float64, model::SemiAnalyticPlummer)::Float64
    x = a / radial_scale(model) 
    return sqrt(1 - 3 / (4 + x^2))
end

function _β_from_α_circular(
    α::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    # Circular orbits: 
    # α = (1 - β^2)^(3/4) / (2 * 3^(3/4) * β^(5/2))
    # use bisection on this function directly
    rootequation(β::Float64) = (1 - β^2)^(3/4) / (2 * 3^(3/4) * β^(5/2)) - α
    return _bisection(rootequation, 1/2, 1)
end