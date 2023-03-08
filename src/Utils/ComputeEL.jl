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
energy as a function of (a,e) for a given potential ψ (and its derivatives)
INCLUDING third derivative of the potential
"""
function EFromAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                 a::Float64,e::Float64,
                 params::OrbitsParameters)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    tole = EccentricityTolerance(a,params.TOLA,params.TOLECC)
    if e <= tole
        # switch to the expanded case
        return EcircExpansion(ψ,dψ,d2ψ,d3ψ,a,e)
    elseif (e == 1.) # up to numerical precision
        return Erad(ψ,a)
    else
        # the analytic version of the energy
        return ((1+e)^(2)*ψ(a*(1+e)) - (1-e)^(2)*ψ(a*(1-e))) / (4e)
    end
end

"""
angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)
"""
function LFromAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                 a::Float64,e::Float64,
                 params::OrbitsParameters)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    tole = EccentricityTolerance(a,params.TOLA,params.TOLECC)
    if e <= tole
        # switch to the expanded case
        return LcircExpansion(dψ,d3ψ,a,e)
    elseif (e == 1.) # up to numerical precision
        return 0.0
    else
        # the analytic version of the angular momentum
        return a * (1-(e)^(2)) * sqrt( (ψ(a*(1+e)) - ψ(a*(1-e))) / (2e) )
    end
end

"""
combined energy + angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)
"""
function ELFromAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                  a::Float64,e::Float64,
                  params::OrbitsParameters)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    E = EFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
    L = LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)

    return E, L
end

########################################################################
#
# (a,e) -> (E,L) mapping : radial values
#
########################################################################

function Erad(ψ::Function,
              a::Float64)::Float64

    return ψ(2*a)
end

########################################################################
#
# (a,e) -> (E,L) mapping : circular Taylor expansions
#
########################################################################

"""
Circular orbit E value
"""
function Ecirc(ψ::F0,dψ::F1,
               a::Float64)::Float64 where {F0 <: Function, F1 <: Function}

    return ψ(a) + 0.5*a*dψ(a)
end

"""
Second-order expansion of energy equation near a circular orbit
"""
function EcircExpansion(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                        a::Float64,e::Float64)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # compute the Taylor expansion of E
    return (0.5*a*dψ(a) + ψ(a)) + (0.5*a*dψ(a) + 0.5*(a)^(2)*d2ψ(a) + (a)^(3)*d3ψ(a)/12) * (e)^(2)
end

"""
Circular orbit L value
"""
function Lcirc(dψ::Function,
               a::Float64)::Float64

    return (sqrt(a))^(3)*sqrt(dψ(a))
end

"""
Coefficients of the second-order expansion of angular momentum equation near a circular orbit
"""
function Lcirc2ndorderExpansionCoefs(dψ::F1,d3ψ::F3,
                                     a::Float64)::Tuple{Float64,Float64,Float64} where {F1 <: Function, F3 <: Function}

    lc = Lcirc(dψ,a)
    coef0 = lc
    coef1 = 0.
    coef2 = ((a)^(5)*d3ψ(a)/(12*lc) - lc)
    return coef0, coef1, coef2
end

"""
Second-order expansion of angular momentum equation near a circular orbit
"""
function LcircExpansion(dψ::F1,d3ψ::F3,
                        a::Float64,e::Float64)::Float64 where {F1 <: Function, F3 <: Function}

    # compute the Taylor expansion of L
    zeroorder, firstorder, secondorder = Lcirc2ndorderExpansionCoefs(dψ,d3ψ,a)
    return zeroorder + firstorder * e + secondorder * (e)^(2)
end


########################################################################
#
# (a,e) -> (E,L) mapping : derivatives, generic case
#
########################################################################

"""
energy and angular momentum derivatives w.r.t. (a,e)
"""
function dELFromAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                   a::Float64,e::Float64,
                   params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)

    tole = EccentricityTolerance(a,params.TOLA,params.TOLECC)
    if e <= tole
        # switch to the expanded case
        ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e = dELcircExpansion(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)

        return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e
    else
        # the analytic version of the energy and angular momentum derivatives w.r.t. (a,e)
        rp, ra = RpRaFromAE(a,e)
        ψrp, ψra, dψrp, dψra = ψ(rp), ψ(ra), dψ(rp), dψ(ra)

        # Difference between potential at apocenter and pericenter
        ψdiff = ψra - ψrp

        dLdenom = 2*sqrt(2*e*ψdiff)

        if (dLdenom == 0.) || isnan(dLdenom) || isinf(dLdenom)
            ∂E∂a, ∂E∂e, ∂L∂a, ∂L∂e = dELcircExpansion(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)
            return E, L, ∂E∂a, ∂E∂e, ∂L∂a, ∂L∂e
        end

        ∂E∂a = ((1+e)^(3)*dψra - (1-e)^(3)*dψrp) / (4e)
        ∂E∂e = (((e)^(2)-1)*ψdiff + a*e*(1+e)^(2)*dψra + a*e*(1-e)^(2)*dψrp) / (4*(e)^(2))

        ∂L∂a = (1-(e)^(2)) * (2*ψdiff + ra*dψra - rp*dψrp)  / (dLdenom)
        ∂L∂e = - a * ((1+3*(e)^(2))*ψdiff - a*e*(1-(e)^(2))*(dψra + dψrp)) / (e*dLdenom)

        return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e
    end
end

"""
energy and angular momentum derivatives w.r.t. (a,e)
EXCLUDING fourth derivative
"""
function dELFromAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                   a::Float64,e::Float64,
                   params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

   # define a numerical fourth derivative
   FDIFF = params.FDIFF
   d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/FDIFF

    E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e = dELFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e
end

"""
energy and angular momentum derivatives w.r.t. (a,e)
EXCLUDING third derivative
"""
function dELFromAE(ψ::Function,dψ::Function,d2ψ::Function,
                   a::Float64,e::Float64,
                   params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

   # define a numerical third derivative
   FDIFF = params.FDIFF
   d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF
   # zero out fourth derivative
   d4ψ(x::Float64) = 0.

    E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e = dELFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e
end


########################################################################
#
# (a,e) -> (E,L) mapping : derivatives, circular Taylor expansions
#
########################################################################

"""
Second-order of energy and angular momentum derivatives w.r.t. (a,e) near circular orbits.
"""
function dELcircExpansion(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                          a::Float64,e::Float64)::Tuple{Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # compute all potential and derivative values
    ψa, dψa, ddψa, dddψa, ddddψa = ψ(a), dψ(a), d2ψ(a), d3ψ(a), d4ψ(a)

    # compute taylor expansion of partial derivatives for E
    ∂E∂a = 0.5*(3*dψa + a*ddψa) + (0.5*dψa + 1.5*a*ddψa + 0.75*(a)^(2)*dddψa + (a)^(3)*ddddψa/12) * (e)^(2)
    ∂E∂e = (a*dψa + (a)^(2)*ddψa + (a)^(3)*dddψa/6) * e

    sqa   = sqrt(a)
    sqdψa = sqrt(dψa)

    # compute taylor expansion of partial derivatives for L
    ∂L∂a = 0.5*(3*sqrt(a)*sqdψa + (sqa)^(3)*ddψa/sqdψa) - 0.5 * sqa * (3*(dψa)^(2) + a*dψa*ddψa - 7*(a)^(2)*dψa*dddψa/12 + (a)^(3)*(ddψa*dddψa - 2*dψa*ddddψa)/12) / ((sqdψa)^(3)) * (e)^(2)
    ∂L∂e = (sqa)^(3) * ((a)^(2)*dddψa/(6*sqdψa) - 2*sqdψa) * e

    return ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e
end


########################################################################
#
# (a,e) -> (E,L) mapping : Jacobian, generic case
#
########################################################################

"""
Jacobian of the (a,e) ↦ (E,L) mapping, i.e. |∂(E,L)/∂(a,e)|
"""
function JacAEToEL(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                   a::Float64,e::Float64,
                   params::OrbitsParameters)::Float64

    _, _, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e = dELFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    return abs(∂E∂a*∂L∂e - ∂L∂a*∂E∂e)
end


"""
Jacobian of the (a,e) ↦ (E,L) mapping, i.e. |∂(E,L)/∂(a,e)|
EXCLUDING fourth derivative
"""
function JacAEToEL(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                   a::Float64,e::Float64,
                   params::OrbitsParameters)::Float64

    FDIFF = params.FDIFF
    d4ψ(r::Float64)::Float64 = (d3ψ(r+FDIFF)-d3ψ(r))/FDIFF

    return JacAEToEL(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end


"""
Jacobian of the (a,e) ↦ (E,L) mapping, i.e. |∂(E,L)/∂(a,e)|
EXCLUDING third derivative (fourth derivative is zer0)
"""
function JacAEToEL(ψ::Function,dψ::Function,d2ψ::Function,
                   a::Float64,e::Float64,
                   params::OrbitsParameters)::Float64

    FDIFF = params.FDIFF
    d3ψ(r::Float64)::Float64 = (d2ψ(r+FDIFF)-d2ψ(r))/FDIFF
    d4ψ(r::Float64)::Float64 = 0.

    return JacAEToEL(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end



########################################################################
#
# (rp,ra) -> (E,L) mapping : from (a,e) -> (E,L) mapping
#
########################################################################

"""
energy as a function of (rp,ra) for a given potential ψ (and its derivatives)
INCLUDING third derivative of the potential
"""
function EFromRpRa(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                   rp::Float64,ra::Float64,
                   params::OrbitsParameters)::Float64

    a,e = AEFromRpRa(rp,ra)

    return EFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
end

"""
energy as a function of (rp,ra) for a given potential ψ (and its derivatives)
EXCLUDING third derivative of the potential
"""
function EFromRpRa(ψ::Function,dψ::Function,d2ψ::Function,
                   rp::Float64,ra::Float64,
                   params::OrbitsParameters)::Float64

    a,e = AEFromRpRa(rp,ra)

    # define a numerical third derivative
    FDIFF = params.FDIFF
    d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF

    return EFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
end

"""
angular momentum as a function of (rp,ra) for a given potential ψ (and its derivatives)
INCLUDING third derivative of the potential
"""
function LFromRpRa(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                   rp::Float64,ra::Float64,
                   params::OrbitsParameters)::Float64

    a,e = AEFromRpRa(rp,ra)

    return LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
end

"""
angular momentum as a function of (rp,ra) for a given potential ψ (and its derivatives)
EXCLUDING third derivative of the potential
"""
function LFromRpRa(ψ::Function,dψ::Function,d2ψ::Function,
                   rp::Float64,ra::Float64,
                   params::OrbitsParameters)::Float64

    a,e = AEFromRpRa(rp,ra)

    # define a numerical third derivative
    FDIFF = params.FDIFF
    d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF

    return LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
end


"""
combined energy + angular momentum as a function of (rp,ra) for a given potenial ψ (and its derivatives)
INCLUDING third derivative
"""
function ELFromRpRa(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                    rp::Float64,ra::Float64,
                    params::OrbitsParameters)::Tuple{Float64,Float64}

    a,e = AEFromRpRa(rp,ra)

    return ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
end

"""
combined energy + angular momentum as a function of (rp,ra) for a given potenial ψ (and its derivatives)
EXCLUDING third derivative
"""
function ELFromRpRa(ψ::Function,dψ::Function,d2ψ::Function,
                    rp::Float64,ra::Float64,
                    params::OrbitsParameters)::Tuple{Float64,Float64}

    a,e = AEFromRpRa(rp,ra)

    # define a numerical third derivative
    FDIFF = params.FDIFF
    d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF

    return ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
end


########################################################################
#
# Circular Radius as a function of angular momentum (mapping inversion)
#
########################################################################


"""RcircFromL(L,dψ,d2ψ[, rmin, rmax])
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
        rcirc = try bisection(r -> L - Lcirc(dψ,r),rmin,rmax,tolx,tolf) catch;   -1. end

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