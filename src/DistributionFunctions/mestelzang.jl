using SpecialFunctions

#####
#
#   The tapered Mestel disc (Zang, see, e.g., Sellwood (2012), Fouvry+ (2015))
#
#####

#####
# Constants dependencies Mestel DF
#####

IntorFloat = Union{Int64,Float64}

"""
    σMestelDF([R0, V0, q])

radial velocity dispersion of the tapered Mestel DF
"""
function σMestelDF(R0::Float64=20.,V0::Float64=1.,q::IntorFloat=11.44)::Float64
    return V0 / sqrt(q+1)
end

"""
    NormConstMestelDF([R0, V0, q])

normalization constant of the tapered Mestel DF.
"""
function NormConstMestelDF(G::Float64,R0::Float64=20.,V0::Float64=1.,q::IntorFloat=11.44)::Float64
    σ = σMestelDF(R0,V0,q)
    return (V0)^(2) / ( 2^(q/2+1) * (pi)^(3/2) * G * gamma(0.5+0.5*q) * (σ)^(q+2) * (R0)^(q+1) )
end

#####
# Initial Mestel DF and derivatives
#####

"""
    MestelDF(E, L[, C, q, sigma])
Mestel distribution function.
"""
function MestelDF(E::Float64,L::Float64,C::Float64=1.444e-14,q::IntorFloat=11.44,σ::Float64=2.835e-1)::Float64

    return C * (L)^(q) * exp(-E / (σ^2))
end
"""
    MesteldDFdE(E, L[, C, q, sigma])
Mestel DF derivative w.r.t. E.
"""
function MesteldDFdE(E::Float64,L::Float64,C::Float64=1.444e-14,q::IntorFloat=11.44,σ::Float64=2.835e-1)::Float64

    return - MestelDF(E,L,C,q,σ) / (σ^2)
end
"""
    MesteldDFdL(E, L[, C, q, sigma])
Mestel DF derivative w.r.t. E.
"""
function MesteldDFdL(E::Float64,L::Float64,C::Float64=1.444e-14,q::IntorFloat=11.44,σ::Float64=2.835e-1)::Float64

    return C * q * (L)^(q-1) * exp(-E / (σ^2))
end

#####
# Zang tapering and derivatives
#####

"""
    ZangInnerTaper(L[, Rin, V0, nu])
Zang inner tapering.
"""
function ZangInnerTaper(L::Float64,Rin::Float64=1.,V0::Float64=1.,ν::Int64=4)::Float64

    return (L^ν) / ( (Rin*V0)^(ν) + (L^ν) )
end
"""
    ZangInnerTaperdL(L[, Rin, V0, nu])
Zang inner tapering derivative.
"""
function ZangInnerTaperdL(L::Float64,Rin::Float64=1.,V0::Float64=1.,ν::Int64=4)::Float64

    return (Rin*V0)^(ν) * ν * (L)^(ν-1) / ( (Rin*V0)^(ν) + (L)^(ν) )^(2)
end

"""
    Zang_outer_tapering(L[, Rout, V0, mu])
Zang outer tapering.
"""
function ZangOuterTaper(L::Float64,Rout::Float64=11.5,V0::Float64=1.,μ::Int64=5)::Float64

    return (Rout*V0)^(μ) / ( (Rout*V0)^(μ) + (L^μ) )
end
"""
    Zang_outer_tapering_dL(L[, Rout, V0, mu])
Zang outer tapering derivative.
"""
function ZangOuterTaperdL(L::Float64,Rout::Float64=11.5,V0::Float64=1.,μ::Int64=5)::Float64

    return - (Rout*V0)^(μ) * μ * (L)^(μ-1) / ( (Rout*V0)^(μ) + (L)^(μ) )^(2)
end

#####
# Full tapered DF and derivatives
#####

"""
    ZangDF(E, L[, R0, Rin, Rout, Rmax, V0, C, q, σ, μ, ν])
Zang star distribution function.
"""
function ZangDF(E::Float64,L::Float64,
                R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                V0::Float64=1.,C::Float64=1.444e-14,
                q::IntorFloat=11.44, σ::Float64=2.835e-1,
                μ::Int64=5, ν::Int64=4)::Float64

    return MestelDF(E,L,C,q,σ) * ZangOuterTaper(L,Rout,V0,μ) * ZangInnerTaper(L,Rin,V0,ν)
end
"""
    ZangdDFdE(E, L[, R0, Rin, Rout, Rmax, V0, C, q, σ, μ, ν])
Zang star DF derivative w.r.t. E.
"""
function ZangdDFdE(E::Float64,L::Float64,
                    R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                    V0::Float64=1.,
                    C::Float64=1.444e-14,
                    q::IntorFloat=11.44, σ::Float64=2.835e-1,
                    μ::Int64=5, ν::Int64=4)::Float64

    return MesteldDFdE(E,L,C,q,σ) * ZangOuterTaper(L,Rout,V0,μ) * ZangInnerTaper(L,Rin,V0,ν)
end
"""
    ZangdDFdL(E, L[, R0, Rin, Rout, Rmax, V0, C, q, σ, μ, ν])
Zang star DF derivative w.r.t. L.
"""
function ZangdDFdL(E::Float64,L::Float64,
                    R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                    V0::Float64=1.,C::Float64=1.444e-14,
                    q::IntorFloat=11.44, σ::Float64=2.835e-1,
                    μ::Int64=5, ν::Int64=4)::Float64

    mesDF = MestelDF(E,L,C,q,σ)
    intap = ZangInnerTaper(L,Rin,V0,ν)
    outap = ZangOuterTaper(L,Rout,V0,μ)
    return ( MesteldDFdL(E,L,C,q,σ) * intap * outap + 
            mesDF * ZangInnerTaperdL(L,Rin,V0,ν) * outap + 
            mesDF * intap * ZangOuterTaperdL(L,Rout,V0,μ) )
end

"""
    ZangndDFdJ(n1, n2, E, L, ndotOmega[, R0, Rin, Rout, Rmax, V0, C, q, σ, μ, ν])
Zang star DF derivative w.r.t. the actions J.
"""
function ZangndDFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,
                    R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                    V0::Float64=1., C::Float64=1.444e-14,
                    q::IntorFloat=11.44, σ::Float64=2.835e-1,
                    μ::Int64=5, ν::Int64=4)::Float64

    if L <= 0.
        println("WARNING: L <= 0.")
        return 0.
    end

    dDFdE = ZangdDFdE(E,L,R0,Rin,Rout,Rmax,V0,C,q,σ,μ,ν)
    dDFdL = ZangdDFdL(E,L,R0,Rin,Rout,Rmax,V0,C,q,σ,μ,ν)
    
    return ndotOmega*dDFdE + n2*dDFdL
end


#####
# Truncated full tapered DF and derivatives
#
# Adding the truncation: no particles beyond Rmax
#####

"""
    TruncatedZangDF(E, L[, R0, Rin, Rout, Rmax, V0, C, q, σ, μ, ν])

Zang star distribution function enforcing ra <= Rmax.
"""
function TruncatedZangDF(E::Float64,L::Float64,
                         R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                         V0::Float64=1.,C::Float64=1.444e-14,
                         q::IntorFloat=11.44, σ::Float64=2.835e-1,
                         μ::Int64=5, ν::Int64=4)::Float64

    if (L <= 0.) || (L > Rmax*V0) || (E < (V0^2)/2 + ψMestel(L/V0,R0,V0)) || (E > ψMestel(Rmax,R0,V0) + L^2/(2*Rmax^2))
        return 0.
    end
    return MestelDF(E,L,C,q,σ) * ZangOuterTaper(L,Rout,V0,μ) * ZangInnerTaper(L,Rin,V0,ν)
end

"""
    TruncatedZangndDFdJ(n1, n2, E, L, ndotOmega[, R0, Rin, Rout, Rmax, V0, C, q, σ, μ, ν])

Zang star DF derivative w.r.t. the actions J enforcing ra <= Rmax.
"""
function TruncatedZangndDFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,
                    R0::Float64=20., Rin::Float64=1., Rout::Float64=11.5, Rmax::Float64=20.,
                    V0::Float64=1., C::Float64=1.444e-14,
                    q::IntorFloat=11.44, σ::Float64=2.835e-1,
                    μ::Int64=5, ν::Int64=4)::Float64

    if (L <= 0.) || (L > Rmax*V0) || (E < (V0^2)/2 + ψMestel(L/V0,R0,V0)) || (E > ψMestel(Rmax,R0,V0) + L^2/(2*Rmax^2))
        return 0.
    end

    dDFdE = ZangdDFdE(E,L,R0,Rin,Rout,Rmax,V0,C,q,σ,μ,ν)
    dDFdL = ZangdDFdL(E,L,R0,Rin,Rout,Rmax,V0,C,q,σ,μ,ν)
    
    return ndotOmega*dDFdE + n2*dDFdL
end