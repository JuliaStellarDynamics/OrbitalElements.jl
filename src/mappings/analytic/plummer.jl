#
#
# This follows K.Tep implementation in CARP
#
#

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
    sp, sa = _spsa_from_rpra(rp, ra, model, params) # extremal anomalies on the orbit
    Ẽ = 1 / sp - (sa^2 - 1) / (sa * sp * (sa + sp)) # dimensionless energy
    L̃ = sqrt(2 * (sp^2 - 1) * (sa^2 - 1) / (sa * sp * (sa + sp))) # dimensionless momentum
    return energy_scale(model) * Ẽ, momentum_scale(model) * L̃ 
end

"""    
for Plummer semi-analytical version, can be found by 1D-bisection instead of 
2D inversion scheme.
"""
function ae_from_EL(
    E::Float64,
    L::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)
    sp, sa = _spsa_from_EL(E, L, model, params)
    rp, ra = _rpra_from_spsa(sp, sa, model, params)
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

@QUESTION: is `w` hard-coded as Hénon's anomaly here ?
"""
function _Θ(
    w::Float64,
    a::Float64,
    e::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)
    rp, ra = rpra_from_ae(a, e) # peri and apocentre
    sp, sa = _spsa_from_rpra(rp, ra, model, params) # extremal anomalies
    
    # Analytical expression of (dr/dw)(1/vr), that is always well-posed
    # see equation (F9) in Tep et al. (2022).
    tep_a = sp * (w + 2) * (w - 1)^2 - sa * (w - 2) * (w + 1)^2
    tep_b = sp * (w^3 - 3w + 6) - sa * (w^3 - 3w - 6)

    return (
        3 
        * sqrt(sa * sp * (sa + sp))
        * tep_a^(3/2)
        / sqrt(sa * sp * tep_a + tep_b)
        / sqrt(4 - w^2)
        / (2 * sqrt(2) * frequency_scale(model))
    )
end

# @IMPROVE:  we do not cure the β computation for near radial orbits using logarithmic 
# sampling near pericentre as in K.Tep's CARP algorithm. For this case the tolerance
# on eccentricity should be increased, at least, to 0.1.

########################################################################
#
# (a,e) ↦ actions mapping
#
########################################################################
# With this anomaly the usual radial action computation is not a good idea.
# Rather use the cured Θ definition than dr/dw
function _radial_action_from_ae(
    a::Float64,
    e::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    # Edge cases
    if a == 0 || e == 0
        return 0.0
    end
    # Handling (a,e)-domain edges through interpolation
    # IMPORTANT : has to be before the generic computation
    fun(atmp::Float64, etmp::Float64) = _radial_action_from_ae(atmp, etmp, model, params)
    res = _interpolate_edges_ae(fun, a, e, params)
    if !isnothing(res)
        return res
    end
    # Generic computations
    function action_integrand(w::Float64)::Float64
        drdw_over_vrad = _Θ(w, a, e, model, params)
        vrad = radial_velocity(w, a, e, model, params)
        return drdw_over_vrad * vrad^2
    end
    return (1 / pi) * _integrate_simpson(action_integrand, params.NINT)
end

########################################################################
#
# (a, e) ↔ frequencies mappings: circular
#
########################################################################
function _Ω1circular(
    r::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    x = r / radial_scale(model) # dimensionless radius
    s = _s_from_r(r, model, params)
    return frequency_scale(model) * sqrt(4 + x^2) / (2 * s^(5/2))
end

function _Ω2circular(
    r::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    s = _s_from_r(r, model)
    return frequency_scale(model) / (2 * s^(3/2))
end

function _βcircular(
    r::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    x = r / radial_scale(model) 
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