using SpecialFunctions

#####
#
#   The tapered Mestel disc (Zang, see, e.g., Sellwood (2012), Fouvry+ (2015))
#
#####

#####
# Constants dependencies Mestel DF
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

#####
# Initial Mestel DF and derivatives
#####

"""
    mestel_DF(E, L[, C, q, sigma])
Mestel distribution function.
"""
function mestel_DF(E::Float64,L::Float64,C::Float64=9.075e-14,q::Float64=11.44,sigma::Float64=2.835e-1)

    return C * (L)^(q) * exp(-E / (sigma)^(2))
end
"""
    mestel_dDFdE(E, L[, C, q, sigma])
Mestel DF derivative w.r.t. E.
"""
function mestel_dDFdE(E::Float64,L::Float64,C::Float64=9.075e-14,q::Float64=11.44,sigma::Float64=2.835e-1)

    return - mestel_DF(E,L,C,q,sigma) / (sigma)^(2)
end
"""
    mestel_dDFdL(E, L[, C, q, sigma])
Mestel DF derivative w.r.t. E.
"""
function mestel_dDFdL(E::Float64,L::Float64,C::Float64=9.075e-14,q::Float64=11.44,sigma::Float64=2.835e-1)

    return C * q * (L)^(q-1.0) * exp(-E / (sigma)^(2))
end

#####
# Zang tapering and derivatives
#####

"""
    Zang_inner_tapering(L[, Rin, V0, nu])
Zang inner tapering.
"""
function Zang_inner_tapering(L::Float64,Rin::Float64=1.,V0::Float64=1.,nu::Int64=4)

    return (L)^(nu) / ( (Rin*V0)^(nu) + (L)^(nu) )
end
"""
    Zang_inner_tapering_dL(L[, Rin, V0, nu])
Zang inner tapering derivative.
"""
function Zang_inner_tapering_dL(L::Float64,Rin::Float64=1.,V0::Float64=1.,nu::Int64=4)

    return (Rin*V0)^(nu) * nu * (L)^(nu-1) / ( (Rin*V0)^(nu) + (L)^(nu) )^(2)
end

"""
    Zang_outer_tapering(L[, Rout, V0, mu])
Zang outer tapering.
"""
function Zang_outer_tapering(L::Float64,Rout::Float64=11.5,V0::Float64=1.,mu::Int64=5)

    return (Rout*V0)^(mu) / ( (Rout*V0)^(mu) + (L)^(mu) )
end
"""
    Zang_outer_tapering_dL(L[, Rout, V0, mu])
Zang outer tapering derivative.
"""
function Zang_outer_tapering_dL(L::Float64,Rout::Float64=11.5,V0::Float64=1.,mu::Int64=5)

    return - (Rout*V0)^(mu) * mu * (L)^(mu-1) / ( (Rout*V0)^(mu) + (L)^(mu) )^(2)
end

#####
# Full tapered DF and derivatives
#####

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
    mestel_Zang_dDFdE(E, L[, R0, Rin, Rout, Rmax, V0, xi, C, q, sigma, nu, mu])
Zang star DF derivative w.r.t. E.
"""
function mestel_Zang_dDFdE(E::Float64,L::Float64;
                        R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                        V0::Float64=1.,
                        xi::Float64=0.5, C::Float64=9.075e-14,
                        q::Float64=11.44, sigma::Float64=2.835e-1,
                        nu::Int64=4, mu::Int64=5)

    return xi * mestel_dDFdE(E,L,C,q,sigma) * Zang_inner_tapering(L,Rin,V0,nu) * Zang_outer_tapering(L,Rout,V0,mu)
end
"""
    mestel_Zang_dDFdL(E, L[, R0, Rin, Rout, Rmax, V0, xi, C, q, sigma, nu, mu])
Zang star DF derivative w.r.t. L.
"""
function mestel_Zang_dDFdL(E::Float64,L::Float64;
                        R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                        V0::Float64=1.,
                        xi::Float64=0.5, C::Float64=9.075e-14,
                        q::Float64=11.44, sigma::Float64=2.835e-1,
                        nu::Int64=4, mu::Int64=5)

    mesDF = mestel_DF(E,L,C,q,sigma)
    intap = Zang_inner_tapering(L,Rin,V0,nu)
    outap = Zang_outer_tapering(L,Rout,V0,mu)
    return xi * (
        mestel_dDFdL(E,L,C,q,sigma) * intap * outap
        + mesDF * Zang_inner_tapering_dL(L,Rin,V0,nu) * outap
        + mesDF * intap * Zang_outer_tapering_dL(L,Rout,V0,mu)
    )
end

"""
    mestel_Zang_ndDFdJ(n1, n2, E, L, ndotOmega[, R0, Rin, Rout, Rmax, V0, xi, C, q, sigma, nu, mu])
Zang star DF derivative w.r.t. the actions J.
"""
function mestel_Zang_ndDFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;
                        R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                        V0::Float64=1.,
                        xi::Float64=0.5, C::Float64=9.075e-14,
                        q::Float64=11.44, sigma::Float64=2.835e-1,
                        nu::Int64=4, mu::Int64=5)

    dDFdE = mestel_Zang_dDFdE(E,L;R0=R0,Rin=Rin,Rout=Rout,Rmax=Rmax,V0=V0,xi=xi,C=C,q=q,sigma=sigma,nu=nu,mu=mu)
    dDFdL = mestel_Zang_dDFdL(E,L;R0=R0,Rin=Rin,Rout=Rout,Rmax=Rmax,V0=V0,xi=xi,C=C,q=q,sigma=sigma,nu=nu,mu=mu)
    
    return ndotOmega*dDFdE + n2*dDFdL
end