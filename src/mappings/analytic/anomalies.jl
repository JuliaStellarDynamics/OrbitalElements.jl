
########################################################################
#
# Appropriate reduced coordinates for Isochrone and Plummer potenials
#
########################################################################
const AnalyticIsochronePlummer = Union{AnalyticIsochrone, SemiAnalyticPlummer}
"""
    _s_from_r(r, model[, params])

effective radius.
"""
function _s_from_r(
    r::Float64,
    model::AnalyticIsochronePlummer,
    params::OrbitalParameters=OrbitalParameters()
)
    x = r / radial_scale(model) # dimensionless radius
    return sqrt(1 + x^2)
end

"""
    _r_from_s(s, model[, params])

radius as a function of effective radius.
"""
function _r_from_s(
    s::Float64,
    model::AnalyticIsochronePlummer,
    params::OrbitalParameters=OrbitalParameters()
)
    # Cure for unexpected effective radius (close to 1, otherwise let the sqrt throw 
    # an error)
    tol = 10*eps(Float64)
    if 1 - tol < s < 1
        return 0.0
    end
    x = sqrt(s^2 - 1) # dimensionless radius
    return radial_scale(model) * x
end

"""
    _r_from_s_derivative(s, model[, params])

derivative of radius w.r.t. effective radius.
"""
function _r_from_s_derivative(
    s::Float64,
    model::AnalyticIsochronePlummer,
    params::OrbitalParameters=OrbitalParameters()
)
    # Cure for unexpected effective radius (close to 1, otherwise let the sqrt throw 
    # an error)
    tol = 10*eps(Float64)
    if 1 - tol < s < 1
        return 0.0
    end
    dxds = s / sqrt(s^2 - 1) # dimensionless radius derivative
    return radial_scale(model) * dxds
end

"""
    _spsa_from_rpra(rp, ra, model[, params])

extremal effective radii of the orbit with pericentre `rp` and apocentre `ra`.
"""
function _spsa_from_rpra(
    rp::Float64,
    ra::Float64,
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)
    # Broadcast the anomaly over peri and apocenter (allocation-free)
    sp, sa = map((r -> _s_from_r(r, model, params)), (rp, ra))
    return sp, sa
end

"""
    _rpra_from_spsa(sp, sa, model::PlummerPotential)

the orbit pericentre `rp` and apocentre `ra` from extremal anomalies.
"""
function _rpra_from_spsa(
    sp::Float64,
    sa::Float64,
    model::SemiAnalyticPlummer,
    params::OrbitalParameters=OrbitalParameters()
)
    # Broadcast the anomaly over peri and apocenter (allocation-free)
    rp, ra = map((s -> _r_from_s(s, model, params)), (sp, sa))
    return rp, ra
end

"""
    _effective_ae_from_spsa(sp, sa)

effective semi major axis and eccentricity from extremal effective radii
"""
function _effective_ae_from_spsa(sp::Float64, sa::Float64)
    # Effective semi-major axis and eccentricity are to sp, sa 
    # what usual semi-major axis and eccentricity are to rp, ra
    return ae_from_rpra(sp, sa)
end

"""
    _effective_ae_from_ae(a, e, model)

effective semi major axis and eccentricity from their usual version
"""
function _effective_ae_from_ae(
    a::Float64,
    e::Float64,
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)
    rp, ra = rpra_from_ae(a, e)
    sp, sa = _spsa_from_rpra(rp, ra, model, params)
    return _effective_ae_from_spsa(sp, sa)
end

function radius_from_anomaly(
    w::Float64,
    a::Float64,
    e::Float64,
    model::AnalyticIsochronePlummer,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    a_eff, e_eff = _effective_ae_from_ae(a, e, model, params)
    s = a_eff * (1 + e_eff * _henonf(w))
    return _r_from_s(s, model, params)
end

function radius_from_anomaly_derivative(
    w::Float64,
    a::Float64,
    e::Float64,
    model::AnalyticIsochronePlummer,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    a_eff, e_eff = _effective_ae_from_ae(a, e, model, params)
    s = a_eff * (1 + e_eff * _henonf(w))
    dsdw = a_eff * e_eff * _henondf(w)
    # Chain rule: dr/dw = dr/ds * ds/dw
    return _r_from_s_derivative(s, model, params) * dsdw
end