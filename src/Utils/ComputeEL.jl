"""Definitions for handling energy and angular momentum

Strategies:
-compute energy and angular momentum from the definitions for (a,e)
-compute derivatives of energy and angular momentum
-add switches for near-circular orbits (E,L,dE,dL)
-include Jacobian to transform between (E,L) and (a,e)
-auxiliary functions to do same calculations as a function of (rp,ra)

"""

########################################################################
#
# Eccentricity tolerance as a function of semi-major axis
#
########################################################################
"""
    EccentricityTolerance(a,TOLA,TOLECC)

tweak eccentricity tolerance to semi-major axis
"""
function EccentricityTolerance(a::Float64,TOLA::Float64,TOLECC::Float64)::Float64

    # We effectively want to switch at
    # a*e > TOLA * TOLECC
    # i.e. tole = TOLA * TOLECC / a
    # + constraint eccentricity tolerance to be in [TOLECC,1/2]
    return min(0.5, max(TOLECC, TOLA*TOLECC/a))
    # return TOLECC
end

########################################################################
#
# (a,e) -> (E,L) mapping : generic case
#
########################################################################

"""
    EFromAE(ψ,dψ,a,e,params)

energy as a function of (a,e) for a given potential ψ (and its derivatives)
"""
function EFromAE(ψ::F0,dψ::F1,
                 a::Float64,e::Float64,
                 params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function}

    tole = EccentricityTolerance(a,params.TOLA,params.TOLECC)
    if (e == 0.)
        return Ecirc(ψ,dψ,a)
    elseif (e == 1.)
        return Erad(ψ,a)
    elseif (0. < e < tole)
        # 2nd order interpolation
        # between circular (e=0) and value at e=tole and e=2*tole
        # Circular:
        ecirc = 0.
        Ec = EFromAE(ψ,dψ,a,ecirc,params)
        # e=tole:
        etole = tole
        Etole = EFromAE(ψ,dψ,a,etole,params)
        # e=2*tole:
        e2tole = 2*tole
        E2tole = EFromAE(ψ,dψ,a,e2tole,params)

        # Interpolation
        return Interpolation2ndOrder(e,ecirc,Ec,etole,Etole,e2tole,E2tole)

    elseif ((1.0-tole) < e < 1.)
        # 2nd order interpolation
        # between radial (e=1.) and value at e=1-tole and e=1-2*tole
        # Radial:
        erad = 1.
        Er = EFromAE(ψ,dψ,a,erad,params)
        # e=1-tole:
        etole = 1.0-tole
        Etole = EFromAE(ψ,dψ,a,etole,params)
        # e=1-2*tole:
        e2tole = 1.0-2*tole
        E2tole = EFromAE(ψ,dψ,a,e2tole,params)

        # Interpolation
        return Interpolation2ndOrder(e,e2tole,E2tole,etole,Etole,erad,Er)
    else
        # the analytic version of the energy
        return ((1+e)^(2)*ψ(a*(1+e)) - (1-e)^(2)*ψ(a*(1-e))) / (4e)
    end
end

"""
    LFromAE(ψ,dψ,a,e,params)

angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)
"""
function LFromAE(ψ::F0,dψ::F1,
                 a::Float64,e::Float64,
                 params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function}

    tole = EccentricityTolerance(a,params.TOLA,params.TOLECC)
    if (e == 0.)
        return Lcirc(dψ,a)
    elseif (e == 1.)
        return 0.
    elseif (0. < e < tole)
        # 2nd order interpolation
        # between circular (e=0) and value at e=tole and e=2*tole
        # Circular:
        ecirc = 0.
        Lc = LFromAE(ψ,dψ,a,ecirc,params)
        # e=tole:
        etole = tole
        Ltole = LFromAE(ψ,dψ,a,etole,params)
        # e=2*tole:
        e2tole = 2*tole
        L2tole = LFromAE(ψ,dψ,a,e2tole,params)

        # Interpolation
        return Interpolation2ndOrder(e,ecirc,Lc,etole,Ltole,e2tole,L2tole)

    elseif ((1.0-tole) < e < 1.)
        # 2nd order interpolation
        # between radial (e=1.) and value at e=1-tole and e=1-2*tole
        # Radial:
        erad = 1.
        Lr = LFromAE(ψ,dψ,a,erad,params)
        # e=1-tole:
        etole = 1.0-tole
        Ltole = LFromAE(ψ,dψ,a,etole,params)
        # e=1-2*tole:
        e2tole = 1.0-2*tole
        L2tole = LFromAE(ψ,dψ,a,e2tole,params)

        # Interpolation
        return Interpolation2ndOrder(e,e2tole,L2tole,etole,Ltole,erad,Lr)
    else
        # the analytic version of the energy
        return a * (1-(e)^(2)) * sqrt( (ψ(a*(1+e)) - ψ(a*(1-e))) / (2e) )
    end
end

"""
    ELFromAE(ψ,dψ,a,e,params)

combined energy + angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)
"""
function ELFromAE(ψ::F0,dψ::F1,
                  a::Float64,e::Float64,
                  params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function}

    E = EFromAE(ψ,dψ,a,e,params)
    L = LFromAE(ψ,dψ,a,e,params)

    return E, L
end

########################################################################
#
# (a,e) -> (E,L) mapping : radial values
#
########################################################################

"""
    Erad(ψ,a)

energy for an exactly radial orbit
"""
function Erad(ψ::Function,
              a::Float64)::Float64

    return ψ(2*a)
end

########################################################################
#
# (a,e) -> (E,L) mapping : circular values
#
########################################################################

"""
    Ecirc(ψ,dψ,a)

energy for an exactly circular orbit.
"""
function Ecirc(ψ::F0,dψ::F1,
               a::Float64)::Float64 where {F0 <: Function, F1 <: Function}

    return ψ(a) + 0.5*a*dψ(a)
end

"""
    Lcirc(dψ,a)

angular momentum for an exactly circular orbit.
"""
function Lcirc(dψ::Function,
               a::Float64)::Float64

    return (sqrt(a))^(3)*sqrt(dψ(a))
end


########################################################################
#
# (a,e) -> (E,L) mapping : derivatives, generic case
#
########################################################################
"""
    dELFromAE(ψ,dψ,d2ψ,a,e,params)

energy and angular momentum derivatives w.r.t. (a,e)
"""
function dELFromAE(ψ::F0,dψ::F1,d2ψ::F2,
                   a::Float64,e::Float64,
                   params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    E, L = ELFromAE(ψ,dψ,a,e,params)

    tole = EccentricityTolerance(a,params.TOLA,params.TOLECC)


    tole = EccentricityTolerance(a,params.TOLA,params.TOLECC)
    if (e == 0.)

        ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e = dELcirc(dψ,d2ψ,a)
        return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e

    elseif (e == 1.)

        ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e = dELrad(ψ,dψ,a)
        return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e

    elseif (0. < e < tole)
        # 2nd order interpolation
        # between circular (e=0) and value at e=tole and e=2*tole
        # Circular:
        ecirc = 0.
        ∂E∂acirc, ∂L∂acirc, ∂E∂ecirc, ∂L∂ecirc = dELcirc(dψ,d2ψ,a)
        # e=tole:
        etole = tole
        _, _, ∂E∂atole, ∂L∂atole, ∂E∂etole, ∂L∂etole = dELFromAE(ψ,dψ,d2ψ,a,etole,params)
        # e=2*tole:
        e2tole = 2*tole
        _, _, ∂E∂a2tole, ∂L∂a2tole, ∂E∂e2tole, ∂L∂e2tole = dELFromAE(ψ,dψ,d2ψ,a,e2tole,params)

        # Interpolation
        ∂E∂a = Interpolation2ndOrder(e,ecirc,∂E∂acirc,etole,∂E∂atole,e2tole,∂E∂a2tole)
        ∂L∂a = Interpolation2ndOrder(e,ecirc,∂L∂acirc,etole,∂L∂atole,e2tole,∂L∂a2tole)
        ∂E∂e = Interpolation2ndOrder(e,ecirc,∂E∂ecirc,etole,∂E∂etole,e2tole,∂E∂e2tole)
        ∂L∂e = Interpolation2ndOrder(e,ecirc,∂L∂ecirc,etole,∂L∂etole,e2tole,∂L∂e2tole)
        
        return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e

    elseif ((1.0-tole) < e < 1.)
        # 2nd order interpolation
        # between radial (e=1.) and value at e=1-tole and e=1-2*tole
        # Radial:
        erad = 1.
        ∂E∂arad, ∂L∂arad, ∂E∂erad, ∂L∂erad = dELrad(ψ,dψ,a)
        # e=1-tole:
        etole = 1.0-tole
        _, _, ∂E∂atole, ∂L∂atole, ∂E∂etole, ∂L∂etole = dELFromAE(ψ,dψ,d2ψ,a,etole,params)
        # e=1-2*tole:
        e2tole = 1.0-2*tole
        _, _, ∂E∂a2tole, ∂L∂a2tole, ∂E∂e2tole, ∂L∂e2tole = dELFromAE(ψ,dψ,d2ψ,a,e2tole,params)

        # Interpolation
        ∂E∂a = Interpolation2ndOrder(e,e2tole,∂E∂a2tole,etole,∂E∂atole,erad,∂E∂arad)
        ∂L∂a = Interpolation2ndOrder(e,e2tole,∂L∂a2tole,etole,∂L∂atole,erad,∂L∂arad)
        ∂E∂e = Interpolation2ndOrder(e,e2tole,∂E∂e2tole,etole,∂E∂etole,erad,∂E∂erad)
        ∂L∂e = Interpolation2ndOrder(e,e2tole,∂L∂e2tole,etole,∂L∂etole,erad,∂L∂erad)

        # Interpolation
        return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e

    else
        # the analytic version of the energy and angular momentum derivatives w.r.t. (a,e)
        rp, ra = RpRaFromAE(a,e)
        ψrp, ψra, dψrp, dψra = ψ(rp), ψ(ra), dψ(rp), dψ(ra)

        # Difference between potential at apocenter and pericenter
        ψdiff = ψra - ψrp

        dLdenom = 2*sqrt(2*e*ψdiff)

        ∂E∂a = ((1+e)^(3)*dψra - (1-e)^(3)*dψrp) / (4e)
        ∂E∂e = (((e)^(2)-1)*ψdiff + a*e*(1+e)^(2)*dψra + a*e*(1-e)^(2)*dψrp) / (4*(e)^(2))

        # apply the brakes if there is a problem!
        if (dLdenom == 0.) || isnan(dLdenom) || isinf(dLdenom)
            return E, L, ∂E∂a, 0., ∂E∂e, 0.
        end

        ∂L∂a = (1-(e)^(2)) * (2*ψdiff + ra*dψra - rp*dψrp)  / (dLdenom)
        ∂L∂e = - a * ((1+3*(e)^(2))*ψdiff - a*e*(1-(e)^(2))*(dψra + dψrp)) / (e*dLdenom)

        return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e
    end
end


########################################################################
#
# (a,e) -> (E,L) mapping : derivatives, circular
#
########################################################################

"""
    dELcirc(ψ,dψ,d2ψ,a)

energy and angular momentum derivatives w.r.t. (a,e) for an exactly circular orbit.
"""
function dELcirc(dψ::F1,d2ψ::F2,
                 a::Float64)::Tuple{Float64,Float64,Float64,Float64} where {F1 <: Function, F2 <: Function}

    dψa, d2ψa = dψ(a), d2ψ(a)
    sqa = sqrt(a)
    
    ∂E∂a = 0.5 * (3.0 * dψa + a*d2ψa)
    ∂E∂e = 0.
    ∂L∂a = 0.5 * (3.0 * sqa * dψa + a * sqa * d2ψa) / sqrt(dψa)
    ∂L∂e = 0.

    return ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e
end


########################################################################
#
# (a,e) -> (E,L) mapping : derivatives, radial
#
########################################################################

"""
    dELrad(ψ,dψ,d2ψ,a)

energy and angular momentum derivatives w.r.t. (a,e) for an exactly radial orbit.
"""
function dELrad(ψ::F0,dψ::F1,
                a::Float64)::Tuple{Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function}

    dψ2a = dψ(2.0*a)

    ∂E∂a = 2.0*dψ2a
    ∂E∂e = a*dψ2a
    ∂L∂a = 0.
    ∂L∂e = - sqrt(2.0) * a * sqrt(ψ(2.0*a)-ψ(0.0))

    return ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e
end

########################################################################
#
# (a,e) -> (E,L) mapping : Jacobian
#
########################################################################

"""
    JacAEToEL(ψ,dψ,d2ψ,a,e,params)

Jacobian of the (a,e) ↦ (E,L) mapping, i.e. |∂(E,L)/∂(a,e)|
"""
function JacAEToEL(ψ::F0,dψ::F1,d2ψ::F2,
                   a::Float64,e::Float64,
                   params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function}

    _, _, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e = dELFromAE(ψ,dψ,d2ψ,a,e,params)

    return abs(∂E∂a*∂L∂e - ∂L∂a*∂E∂e)
end


########################################################################
#
# (rp,ra) -> (E,L) mapping : from (a,e) -> (E,L) mapping
#
########################################################################

"""
    EFromRpRa(ψ,dψ,rp,ra,params)

energy as a function of (rp,ra) for a given potential ψ (and its derivatives)
"""
function EFromRpRa(ψ::F0,dψ::F1,
                   rp::Float64,ra::Float64,
                   params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function}

    a,e = AEFromRpRa(rp,ra)

    return EFromAE(ψ,dψ,a,e,params)
end


"""
    LFromRpRa(ψ,dψ,rp,ra,params)

angular momentum as a function of (rp,ra) for a given potential ψ (and its derivatives)
"""
function LFromRpRa(ψ::F0,dψ::F1,
                   rp::Float64,ra::Float64,
                   params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function}

    a,e = AEFromRpRa(rp,ra)

    return LFromAE(ψ,dψ,a,e,params)
end


"""
    ELFromRpRa(ψ,dψ,rp,ra,params)

combined energy + angular momentum as a function of (rp,ra) for a given potenial ψ (and its derivatives)
"""
function ELFromRpRa(ψ::F0,dψ::F1,
                    rp::Float64,ra::Float64,
                    params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function}

    a,e = AEFromRpRa(rp,ra)

    return ELFromAE(ψ,dψ,a,e,params)
end


########################################################################
#
# Circular Radius as a function of angular momentum (mapping inversion)
#
########################################################################


"""
    RcircFromL(L,dψ,rmin,rmax,tolx,tolf)
    
perform backwards mapping from L for a circular orbit to radius
can tune [rmin,rmax] for extra optimisation (but not needed)
@WARNING: important assumption Ω1circular is a decreasing function of radius
"""
function RcircFromL(L::Float64,
                    dψ::F1,
                    rmin::Float64,rmax::Float64,
                    tolx::Float64=1000.0*eps(Float64),tolf::Float64=1000.0*eps(Float64))::Float64 where {F1 <: Function}

    # check that the input frequency is valid
    if L  < 0.
        error("OrbitalElements.Utils.RcircFromL: Negative angular momentum L = $L")
        return -1.
    elseif L == 0.
        return 0.
    else 
        # use bisection to find the circular orbit radius corresponding to given frequency
        rcirc = try bisection(r -> L - Lcirc(dψ,r),rmin,rmax,tolx=tolx,tolf=tolf) catch;   -1. end

        # check if bisection failed: report why
        if (rcirc == -1.)
            if (Lcirc(dψ,rmax) < L)
                return RcircFromL(L,dψ,rmax,10*rmax,tolx,tolf)
            elseif (L < Lcirc(dψ,rmin))
                return RcircFromL(L,dψ,rmin/10,rmin,tolx,tolf)
            else
                error("OrbitalElements.Utils.RcircFromL: Unable to find the associated radius of L = $L")
                return -1.
            end
        end

        return rcirc
    end
end