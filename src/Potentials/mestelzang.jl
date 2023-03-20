#####
#
#   The truncated Mestel disc (Zang, see, e.g., Sellwood (2012), Fouvry+ (2015))
#   + potential flattening in the center
#
#####

#####
#
# Classic Mestel potential
#
#####

"""
    ψMestel(r[, R0, V0, epsilon])

the Mestel potential (flat rotation curve).
"""
function ψMestel(r::Float64,R0::Float64=1.,V0::Float64=1.)

    return (V0)^(2) * log(r/R0)
end

"""
    dψMestel(r[, R0, V0, epsilon])

the Mestel potential derivative.
"""
function dψMestel(r::Float64,R0::Float64=1.,V0::Float64=1.)

    x = r/R0
    return ((V0)^(2) / R0) / x
end

"""
    d2ψMestel(r[, R0, V0, epsilon])

the Mestel potential second derivative.
"""
function d2ψMestel(r::Float64,R0::Float64=1.,V0::Float64=1.)

    x = r/R0
    return  - ((V0)^(2) / (R0)^(2)) / (x^2)
end
"""
    d3ψMestel(r[, R0, V0, epsilon])

the Mestel potential third derivative.
"""
function d3ψMestel(r::Float64,R0::Float64=1.,V0::Float64=1.)

    x = r/R0
    return  2.0 * ((V0)^(2) / (R0)^(3)) / (x^3)
end
"""
    d4ψMestel(r[, R0, V0, epsilon])

the Mestel potential fourth derivative.
"""
function d4ψMestel(r::Float64,R0::Float64=1.,V0::Float64=1.)

    x = r/R0
    return  - 6.0 * ((V0)^(2) / (R0)^(4)) / (x^4)
end

"""
    Ω₀Mestel([R0, V0, epsilon])

the truncated Mestel frequency scale.
"""
function Ω₀Mestel(R0::Float64=1.,V0::Float64=1.)

    return 2.0 * V0 / R0
end

#####
#
# Tapered Mestel potential
#
#####

"""
    ψMestelTrunc(r[, R0, V0, epsilon])

the truncated Mestel potential (flat rotation curve).
"""
function ψMestelTrunc(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    return 0.5 * (V0)^(2) * log((r/R0)^(2) + (eps)^(2))
end

"""
    dψMestelTrunc(r[, R0, V0, epsilon])

the truncated Mestel potential derivative.
"""
function dψMestelTrunc(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    x = r/R0
    return ((V0)^(2) / R0) * x / ((eps)^(2) + (x)^(2))
end

"""
    d2ψMestelTrunc(r[, R0, V0, epsilon])

the truncated Mestel potential second derivative.
"""
function d2ψMestelTrunc(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    x = r/R0
    return  ((V0)^(2) / (R0)^(2)) * ((eps)^(2) - (x)^(2)) / ((eps)^(2) + (x)^(2))^(2)
end
"""
    d3ψMestelTrunc(r[, R0, V0, epsilon])

the truncated Mestel potential third derivative.
"""
function d3ψMestelTrunc(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    x = r/R0
    return  ((V0)^(2) / (R0)^(3)) * 2 * x * ((x)^2 - 3*(eps)^(2)) / ((eps)^(2) + (x)^(2))^(3)
end
"""
    d4ψMestelTrunc(r[, R0, V0, epsilon])

the truncated Mestel potential fourth derivative.
"""
function d4ψMestelTrunc(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    x = r/R0
    return  - 6 * ((V0)^(2) / (R0)^(4)) * ((x)^(4) - 6*(x)^(2)*(eps)^(2) + (eps)^(4)) / ((eps)^(2) + (x)^(2))^(4)
end

"""
    Ω₀MestelTrunc([R0, V0, epsilon])

the truncated Mestel frequency scale.
"""
function Ω₀MestelTrunc(R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    return (2.0 * V0) / (eps * R0)
end
