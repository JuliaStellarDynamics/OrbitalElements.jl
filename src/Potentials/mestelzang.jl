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

    return (V0)^(2)*log( sqrt( (r/R0)^(2)+(eps)^(2) ) )
end

"""
    mestel_dpsi_dr(r[, R0, V0, epsilon])

the truncated Mestel potential derivative.
"""
function mestel_dpsi_dr(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    return (V0)^(2) / R0 * (r/R0) / ( (eps)^(2) + (r/R0)^(2) )
end

"""
    mestel_ddpsi_ddr(r[, R0, V0, epsilon])

the truncated Mestel potential second derivative.
"""
function mestel_ddpsi_ddr(r::Float64,R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    return - (V0)^(2) / (R0)^(2) * ( (eps)^(2) - (r/R0)^2 ) / ( (eps)^(2) + (r/R0)^(2) )^(2)
end

"""
    mestel_Omega0([R0, V0, epsilon])

the truncated Mestel frequency scale.
"""
function mestel_Omega0(R0::Float64=1.,V0::Float64=1.,eps::Float64=0.01)

    return (2.0 * V0) / (eps * R0)
end


#####
#
# Tapered Mestel DF
#
#####


"""
    sigmar_Mestel_DF([R0, V0, q])

radial velocity dispersion of the tapered Mestel DF
"""
function sigmar_Mestel_DF(R0::Float64=20.,V0::Float64=1.,q::Float64=11.44)
    return V0 / sqrt(q+1.0)
end

"""
    normC_Mestel_DF([R0, V0, q])

normalization constant of the tapered Mestel DF.
"""
function normC_Mestel_DF(R0::Float64=20.,V0::Float64=1.,q::Float64=11.44)
    sigma = sigmar_Mestel_DF(V0,q)
    return (V0)^(2) / ( (2.0)^(0.5*q) * sqrt(pi) * gamma(0.5+0.5*q) * (sigma)^(q+2.0) * (R0)^(q+1.0) )
end


"""
    mestel_DF(E, L[, C, q, sigma])
Mestel distribution function.
"""
function mestel_DF(E::Float64,L::Float64,C::Float64=9.075e-14,q::Float64=11.44,sigma::Float64=2.835e-1)

    return C * (L)^(q) * exp(-E / (sigma)^(2))
end

"""
    Zang_inner_tapering(L[, Rin, V0, nu])
Zang inner tapering.
"""
function Zang_inner_tapering(L::Float64,Rin::Float64=1.,V0::Float64=1.,nu::Int64=4)

    return (L)^(nu) / ( (Rin*V0)^(nu) + (L)^(nu) )
end

"""
    Zang_outer_tapering(L[, Rout, V0, mu])
Zang outer tapering.
"""
function Zang_outer_tapering(L::Float64,Rout::Float64=11.5,V0::Float64=1.,mu::Int64=5)

    return (Rout*V0)^(mu) / ( (Rout*V0)^(mu) + (L)^(mu) )
end

"""
    mestel_Zang_DF(E, L[, R0, Rin, Rout, Rmax, V0, xi, C, q, sigma, nu, mu])
Zang star distribution function.

@ WARNING : Cut off at large radius to add !
"""
function mestel_Zang_DF(E::Float64,L::Float64;
                        R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                        V0::Float64=1.,
                        xi::Float64=0.5, C::Float64=9.075e-14, 
                        q::Float64=11.44, sigma::Float64=2.835e-1,
                        nu::Int64=4, mu::Int64=5)

    return xi * mestel_DF(E,L,C,q,sigma) * Zang_inner_tapering(L,Rin,V0,nu) * Zang_outer_tapering(L,Rout,V0,mu)
end


"""
    mestel_Zang_DF(E, L[, R0, Rin, Rout, Rmax, V0, xi, C, q, sigma, nu, mu])
Zang star distribution function.

@ WARNING : Cut off at large radius to add !
"""
function mestel_Zang_DF(E::Float64,L::Float64;
                        R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                        V0::Float64=1.,
                        xi::Float64=0.5, C::Float64=9.075e-14, 
                        q::Float64=11.44, sigma::Float64=2.835e-1,
                        nu::Int64=4, mu::Int64=5)

    return xi * mestel_DF(E,L,C,q,sigma) * Zang_inner_tapering(L,Rin,V0,nu) * Zang_outer_tapering(L,Rout,V0,mu)
end