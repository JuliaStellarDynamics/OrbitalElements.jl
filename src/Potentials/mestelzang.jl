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
    Ω₀Mestel([R0, V0, epsilon])

the truncated Mestel frequency scale.
"""
function Ω₀Mestel(R0::Float64=1.,V0::Float64=1.)

    return V0 / R0
end

#####
#
# Tapered Mestel potential
#
#####

"""
    ψMestelTrunc(r[, R0, V0, ε0])

the truncated Mestel potential (flat rotation curve).
"""
function ψMestelTrunc(r::Float64,R0::Float64=1.,V0::Float64=1.,ε0::Float64=0.01)

    return 0.5 * (V0)^(2) * log((r/R0)^(2) + (ε0)^(2))
end

"""
    dψMestelTrunc(r[, R0, V0, ε0])

the truncated Mestel potential derivative.
"""
function dψMestelTrunc(r::Float64,R0::Float64=1.,V0::Float64=1.,ε0::Float64=0.01)

    x = r/R0
    # Stable version at infinity (not stable in 0.)
    if x > 1.e5
        return ((V0)^(2) / R0) / (x * (1.0 + (ε0/x)^(2)))
    end

    return ((V0)^(2) / R0) * x / ((ε0)^(2) + (x)^(2))
end

"""
    d2ψMestelTrunc(r[, R0, V0, ε0])

the truncated Mestel potential second derivative.
"""
function d2ψMestelTrunc(r::Float64,R0::Float64=1.,V0::Float64=1.,ε0::Float64=0.01)

    x = r/R0
    # Stable version at infinity (not stable in 0.)
    if x > 1.e5
        return ((V0)^(2) / (R0)^(2)) * ((ε0/x)^(2) - 1.0) / (x * (1.0 + (ε0/x)^(2)))^2
    end

    return  ((V0)^(2) / (R0)^(2)) * ((ε0)^(2) - (x)^(2)) / ((ε0)^(2) + (x)^(2))^(2)
end

"""
    Ω₀MestelTrunc([R0, V0, epsilon])

the truncated Mestel frequency scale.
"""
function Ω₀MestelTrunc(R0::Float64=1.,V0::Float64=1.,ε0::Float64=0.01)

    return V0 / R0
end
