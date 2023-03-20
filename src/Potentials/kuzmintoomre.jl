#####
#
#   The Kuzmin-Toomre disc (Same as Plummer ...)
#
#####

"""ψKT(r[, bc, M, G])

the Kuzmin-Toomre disc potential
"""
function ψKT(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    return -astronomicalG*M*(sqrt(r^2 + bc^2))^(-1)
end

"""dψKT(r[, bc, M, G])

the Kuzmin-Toomre disc potential derivative
"""
function dψKT(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    return astronomicalG*M*r*(sqrt(r^2 + bc^2))^(-3)
end

"""d2ψKT(r[, bc, M, G])

the Kuzmin-Toomre disc potential second derivative
"""
function d2ψKT(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    return astronomicalG*M*(bc^2 - 2*(r^2))*(sqrt(r^2 + bc^2))^(-5)
end

"""
    Ω₀KT([bc, M, G])

the Kuzmin-Toomre disc frequency scale.
"""
function Ω₀KT(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    return 2*sqrt(astronomicalG*M/bc^3)
end
