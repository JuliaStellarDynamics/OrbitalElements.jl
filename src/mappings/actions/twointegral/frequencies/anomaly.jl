"""

@IMPROVE: Henon anomaly is hard-coded. Make it a function of parameters ?
Mainly in 3 places: `_Θ_edge`, `radius_from_anomaly` and `radius_from_anomaly_derivative`
"""

########################################################################
#
# Hénon anomaly function and derivatives
#
########################################################################
"""
    henonf(x)

the henon anomaly increment
"""
_henonf(x::Float64)::Float64 = x * (1.5 - 0.5 * x^2)

"""
    henondf(x)

the derivative of the henon anomaly increment
"""
_henondf(x::Float64)::Float64 = 1.5 * (1.0 - x^2)

"""
    henond2f(x)

the second derivative of the henon anomaly increment
"""
_henond2f(x::Float64)::Float64 = -3x

"""
    henond3f(x)

the third derivative of the henon anomaly increment
"""
_henond3f(x::Float64)::Float64 = -3.

"""
    henond4f(x)

the fourth derivative of the henon anomaly increment
"""
_henond4f(x::Float64)::Float64 = 0.


########################################################################
#
# Radius as a function of u and derivatives in Hénon anomaly mapping
#
########################################################################
"""
    radius_from_anomaly(w, a, e, model[, params])

mapping from anomaly to radius. Default is Henon anomaly.

@IMPROVE: right now, Hénon anomaly is hard-coded. (except for AnalyticIsochrone and 
SemiAnalyticPlummer, for which it is redefined using specific anomaly)
"""
function radius_from_anomaly(
    w::Float64,
    a::Float64,
    e::Float64,
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    return a * (1 + e * _henonf(w))
end

"""
    radius_from_anomaly_derivative(w, a, e, model[, params])

derivative of the mapping from anomaly to radius. Default is Henon anomaly.

@IMPROVE: right now, Hénon anomaly is hard-coded. (except for AnalyticIsochrone and 
SemiAnalyticPlummer, for which it is redefined using specific anomaly)
"""
function radius_from_anomaly_derivative(
    w::Float64,
    a::Float64,
    e::Float64,
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    return a * e * _henondf(w)
end