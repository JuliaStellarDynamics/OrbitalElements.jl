#####
#
#   The truncated Mestel disc (Zang, see, e.g., Sellwood (2012), Fouvry+ (2015))
#   + potential flattening in the center
#
#####


#####
#
# Tapered Mestel potential
#
#####

"""
    mestel_psi(r[, R0, V0, epsilon])

the truncated Mestel potential (flat rotation curve).
"""
function mestel_psi(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    return 0.5 * (V0)^(2) * log((r/R0)^(2) + (eps)^(2))
end

"""
    mestel_dpsi_dr(r[, R0, V0, epsilon])

the truncated Mestel potential derivative.
"""
function mestel_dpsi_dr(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    x = r/R0
    return ((V0)^(2) / R0) * x / ((eps)^(2) + (x)^(2))
end

"""
    mestel_d2psi_dr2(r[, R0, V0, epsilon])

the truncated Mestel potential second derivative.
"""
function mestel_ddpsi_ddr(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    x = r/R0
    return  ((V0)^(2) / (R0)^(2)) * ((eps)^(2) - (x)^(2)) / ((eps)^(2) + (x)^(2))^(2)
end
"""
    mestel_d3psi_dr3(r[, R0, V0, epsilon])

the truncated Mestel potential third derivative.
"""
function mestel_d3psi_dr3(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    x = r/R0
    return  ((V0)^(2) / (R0)^(3)) * 2 * x * ((x)^2 - 3*(eps)^(2)) / ((eps)^(2) + (x)^(2))^(3)
end
"""
    mestel_d4psi_dr4(r[, R0, V0, epsilon])

the truncated Mestel potential fourth derivative.
"""
function mestel_d4psi_dr4(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    x = r/R0
    return  - 6 * ((V0)^(2) / (R0)^(4)) * ((x)^(4) - 6*(x)^(2)*(eps)^(2) + (eps)^(4)) / ((eps)^(2) + (x)^(2))^(4)
end

"""
    mestel_Omega0([R0, V0, epsilon])

the truncated Mestel frequency scale.
"""
function mestel_Omega0(R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    return (2.0 * V0) / (eps * R0)
end