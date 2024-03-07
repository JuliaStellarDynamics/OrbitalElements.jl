"""

@IMPROVE: Henon anomaly is hard-coded. Make it a function of parameters ?
Mainly in 3 places: `Θ_edge`, `ru` and `drdu`
"""

########################################################################
#
# Hénon anomaly function and derivatives
#
########################################################################
"""
    henonf(u)

the henon anomaly increment
"""
_henonf(u::Float64)::Float64 = u * (1.5 - 0.5 * u^2)

"""
    henondf(u)

the derivative of the henon anomaly increment
"""
_henondf(u::Float64)::Float64 = 1.5 * (1.0 - u^2)

"""
    henond2f(u)

the second derivative of the henon anomaly increment
"""
_henond2f(u::Float64)::Float64 = -3u

"""
    henond3f(u)

the third derivative of the henon anomaly increment
"""
_henond3f(u::Float64)::Float64 = -3.

"""
    henond4f(u)

the fourth derivative of the henon anomaly increment
"""
_henond4f(u::Float64)::Float64 = 0.


########################################################################
#
# Radius as a function of u and derivatives in Hénon anomaly mapping
#
########################################################################
"""
    radius_from_anomaly(u, a, e)

mapping from u->r in Henon variables.

@IMPROVE: right now, Hénon anomaly is hard-coded.
"""
radius_from_anomaly(u::Float64,a::Float64,e::Float64)::Float64 = a * (1 + e * _henonf(u))

"""
    radius_from_anomaly_derivative(u, a, e)

mapping from u->r in Henon variables

@IMPROVE: right now, Hénon anomaly is hard-coded.
"""
function radius_from_anomaly_derivative(u::Float64,a::Float64,e::Float64)::Float64
    return a * e * _henondf(u)
end