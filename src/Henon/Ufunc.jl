"""
Ufunc.jl

collection of functions that are used to compute quantities when integrating over an orbit angle.

Specific to Henon mapping.

"""

"""ruRpRa(u,rp,ra)
mapping from u->r in Henon variables
"""
function ruRpRa(u::Float64,rp::Float64,ra::Float64)
    a,e = ae_from_rpra(rp,ra)
    fu = u*(3/2 - (u^2)/2)
    return a*(1+e*fu)
end

"""ru(u,a,e)
mapping from u->r in Henon variables

one might also consider making a version that takes a passable function (i.e. if not using Henon)
function ru(fun::Function,u::Float64,a::Float64,e::Float64)
    return a*(1+e*fun(u))
end

"""
function ru(u::Float64,a::Float64,e::Float64)
    return a*(1+e*henon_f(u))
end

function ruAE(u::Float64,a::Float64,e::Float64)
    return a*(1+e*henon_f(u))
end

"""drduRpRa(u,rp,ra)
first derivative of r(u)
"""
function drduRpRa(u::Float64,rp::Float64,ra::Float64)
    a,e = ae_from_rpra(rp,ra)
    return (3/2)*a*e*(1-(u^2))
end

"""d2rduRpRa(u,rp,ra)
second derivative of r(u)
"""
function d2rduRpRa(u::Float64,rp::Float64,ra::Float64)
    a,e = ae_from_rpra(rp,ra)
    return -3*a*e*u
end

"""d3rduRpRa(u,rp,ra)
third derivative of r(u)
"""
function d3rduRpRa(u::Float64,rp::Float64,ra::Float64)
    a,e = ae_from_rpra(rp,ra)
    return -3*a*e
end

"""
all additional derivatives of r(u) are 0
"""


"""
the effective potential: note the relationship to Q
"""
function ψeff(ψ::Function,r::Float64,L::Float64)::Float64
    return ψ(r) + (1/2) * (L/r)^(2)
end

"""
the derivative of the effective potential
"""
function dψeffdr(dψ::Function,r::Float64,L::Float64)::Float64
    return dψ(r) - (L)^(2) / (r^3)
end

"""
the second derivative of the effective potential
"""
function d2ψeffdr2(d2ψ::Function,r::Float64,L::Float64)::Float64
    return d2ψ(r) + 3 * (L)^(2) / (r^4)
end



"""Q(ψ,dψ/dr,d²ψ/dr²,u,rp,ra)
Q \equiv 2(E-psi(r)) - L^2/r^2
uses mapping r(u)
"""
function Q(ψ::Function,
           dψ::Function,
           d2ψ::Function,
           u::Float64,
           rp::Float64,
           ra::Float64)::Float64

    E = EFromRpRa(ψ,dψ,d2ψ,rp,ra)
    L = LFromRpRa(ψ,dψ,d2ψ,rp,ra)

    r  = ruRpRa(u,rp,ra)

    # if radial orbits, ignore the contribution from the second term
    if L==0
        return 2*(E-ψ(r))
    else
        return 2*(E-ψ(r)) - (L^2)/(r^2)
    end
end

"""ExpandThetaRpRa(ψ,dψ/dr,d²ψ/dr²,u,rp,ra[,VERBOSE])
expansion of Theta near u=pm1.
"""
function ExpandThetaRpRa(ψ::Function,
                         dψ::Function,
                         d2ψ::Function,
                         u::Float64,
                         rp::Float64,
                         ra::Float64;
                         VERBOSE::Int64=0,
                         SECONDORDER::Bool=true)

    # check for radial orbits
    L = LFromRpRa(ψ,dψ,d2ψ,rp,ra)

    if (L==0)
        return 0.
    end

    # check the sign of u
    if u>0
        usign = 1.0
    else
        usign = -1.0
    end

    # derivatives of Q
    dQdrval = dQdr(ψ,dψ,d2ψ,usign,rp,ra)
    ddQddrval = ddQddr(ψ,dψ,d2ψ,usign,rp,ra)

    # derivative of the anomaly
    ddr = d2rduRpRa(usign,rp,ra)
    dddr = d3rduRpRa(usign,rp,ra)

    # if the sqrt is about to be zero...don't let it.
    #if ((ddr < 0) ⊻ (dQdrval < 0)) # XOR gate method, is this faster or slower?
    if (ddr*dQdrval < 0)

        if VERBOSE>0
            println("OrbitalElements/Ufunc.jl: Bad Theta expansion=",ddr,",",dQdrval," (u=",u,",rp=",rp,",ra=",ra,")")
        end

        u = -1.
        ddr = d2rduRpRa(u,rp,ra)
        dQdrval = dQdr(ψ,dψ,d2ψ,u,rp,ra)

    end

    zeroorder   = sqrt(2*ddr/dQdrval)
    firstorder  = -1usign*sqrt(2)*dddr*(u-usign)/(3*sqrt(dQdrval*ddr))
    secondorder = usign*((3ddQddrval*(ddr^3) + dQdrval*(dddr^2))*((u-usign)^2))/(12*sqrt(2)*((dQdrval*ddr)^(3/2)))

    if SECONDORDER
        return zeroorder+firstorder+secondorder
    else
        return zeroorder+firstorder
    end

end


"""Theta(ψ,dψ/dr,d²ψ/dr²,u,rp,ra[,EDGE,VERBOSE])

``Theta(s) \\equiv (1/vr)*(dr/ds) where s\\in[-1,1]``

EDGE sets the switch to the expanded approximation near the boundaries
for radial orbits, the pericentre limit is overridden and set to 0.
"""
function ThetaRpRa(ψ::Function,
                   dψ::Function,
                   d2ψ::Function,
                   u::Float64,
                   rp::Float64,
                   ra::Float64;
                   EDGE::Float64=0.01,
                   VERBOSE::Int64=0,
                   SECONDORDER::Bool=true)

    L = LFromRpRa(ψ,dψ,d2ψ,rp,ra)

    if (L==0)
        EDGE=1.e-10
    end

    if ((1-abs(u))<EDGE)
        # if near the edge, use the expanded approximations
        res = ExpandThetaRpRa(ψ,dψ,d2ψ,u,rp,ra)

    else # main branch, direct calculation

        Qval = Q(ψ,dψ,d2ψ,u,rp,ra)

        # check if Q is negative. If negative...revert back to the expanded values
        if (Qval<0)
            if VERBOSE>0
                println("OrbitalElements/Ufunc.jl: Bad Q=",Qval," (u=",u,",rp=",rp,",ra=",ra,")")
            end

            res = ExpandThetaRpRa(ψ,dψ,d2ψ,u,rp,ra)

        # if not negative, proceed along main branch
        else
            dr  = drduRpRa(u,rp,ra)
            res = (1/sqrt(Qval))*dr
        end

    end

    return res
end

"""dQdr(ψ,dψ/dr,d²ψ/dr²,u,rp,ra)

Partial derivative of Q w.r.t. r

@IMPROVE, watch for switch when L==0

"""
function dQdr(ψ::Function,
              dψ::Function,
              d2ψ::Function,
              u::Float64,
              rp::Float64,
              ra::Float64)

    r = ruRpRa(u,rp,ra)
    L = LFromRpRa(ψ,dψ,d2ψ,rp,ra)

    if (L==0)
        return max(0.,-2*dψ(r)-2*d2ψ(r)*r)
    else
        return  ((2*(L^2)/(r^3)) - 2*dψ(r))
    end

end

"""ddQddr(ψ,dψ/dr,d²ψ/dr²,u,rp,ra)

Second partial derivative of Q w.r.t. r

@IMPROVE, watch for switch when L==0

"""
function ddQddr(ψ::Function,
                dψ::Function,
                d2ψ::Function,
                u::Float64,
                rp::Float64,
                ra::Float64)

    r   = ruRpRa(u,rp,ra)
    L = LFromRpRa(ψ,dψ,d2ψ,rp,ra)

    if (L==0)
        return max(0.,-2*dψ(r)-2*d2ψ(r)*r)
    else
        return  ((-6*(L^2)/(r^4)) - 2*d2ψ(r))
    end

end


"""dThetadu(ψ,dψ/dr,d²ψ/dr²,u,rp,ra)
partial derivative of Theta w.r.t. u.
"""
function dThetadu(ψ::Function,
                  dψ::Function,
                  d2ψ::Function,
                  u::Float64,
                  rp::Float64,
                  ra::Float64)

    Qval = Q(ψ,dψ,d2ψ,u,rp,ra)

    r   = ruRpRa(u,rp,ra)
    dr  = drduRpRa(u,rp,ra)
    ddr = d2rduRpRa(u,rp,ra)

    dQdu = dQdr(ψ,dψ,d2ψ,u,rp,ra)*drduRpRa(u,rp,ra)

    return (Qval^(-1/2))*ddr + dr*((-1/2)*dQdu*(Qval^(-3/2)))

end

"""dQdrnumerical(ψ,dψ/dr,d²ψ/dr²,u,rp,ra[,eps])
numerical check version of dQ/dr
"""
function dQdrnumerical(ψ::Function,
                       dψ::Function,
                       d2ψ::Function,
                       u::Float64,
                       rp::Float64,
                       ra::Float64,
                       eps::Float64=0.001)

    Qval1 = Q(ψ,dψ,d2ψ,u,rp,ra)
    Qval2 = Q(ψ,dψ,d2ψ,u+eps,rp,ra)
    r1   = ruRpRa(u,rp,ra)
    r2   = ruRpRa(u+eps,rp,ra)
    return (Qval2-Qval1)/(r2-r1)
end

"""dQdunumerical(ψ,dψ/dr,d²ψ/dr²,u,rp,ra[,eps])
numerical check version of dQ/du
"""
function dQdunumerical(ψ::Function,
                       dψ::Function,
                       d2ψ::Function,
                       u::Float64,
                       rp::Float64,
                       ra::Float64;
                       eps::Float64=0.001)

    Qval1 = Q(ψ,dψ,d2ψ,u,rp,ra)
    Qval2 = Q(ψ,dψ,d2ψ,u+eps,rp,ra)
    return (Qval2-Qval1)/(eps)

end

"""dThetadunumerical(ψ,dψ/dr,d²ψ/dr²,u,rp,ra)
numerical check version of dTheta/du
"""
function dThetadunumerical(ψ::Function,
                           dψ::Function,
                           d2ψ::Function,
                           u::Float64,
                           rp::Float64,
                           ra::Float64;
                           epse::Float64=0.001)

    Tval1 = Theta(ψ,dψ,d2ψ,u,rp,ra,0.)
    Tval2 = Theta(ψ,dψ,d2ψ,u+eps,rp,ra,0.)
    return (Tval2-Tval1)/(eps)

end



"""

ThetaExpansion, the anomaly for computing orbit averages
as a function of (a,e)

@IMPROVE, once happy with this version, remove ThetaRpRa

"""
function ThetaExpansionAE(ψ::Function,
                          dψ::Function,
                          d2ψ::Function,
                          d3ψ::Function,
                          u::Float64,
                          a::Float64,
                          e::Float64;
                          TOLECC::Float64=ELTOLECC,
                          f::Function=henon_f,
                          df::Function=henon_df,
                          d2f::Function=henon_d2f,
                          d3f::Function=henon_d3f,
                          d4f::Function=henon_d4f)::Float64

    ul = (u > 0) ? 1.0 : -1.0
    #rl = ru(f,ul,a,e)
    rl = ru(ul,a,e)

    E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)

    dψeffl, d2ψeffl = dψeffdr(dψ,rl,L), d2ψeffdr2(d2ψ,rl,L)

    d2fl, d3fl, d4fl = d2f(ul), d3f(ul), d4f(ul)

    pref = - ul * a * e

    # this denominator can be negative?
    combination = - a * e * dψeffl * d2fl

    # switch to safety: don't contribute anything at this point
    if combination <= 0.
        return 0.0
    end

    denom = sqrt(combination)

    zeroorder   = d2fl
    firstorder  = d3fl / 3.0
    secondorder = (3.0 * (dψeffl * d2fl * d4fl - a * e * d2ψeffl * (d2fl)^(3)) -  dψeffl * (d3fl)^(2)) / (24 * dψeffl * d2fl)

    return pref / denom * ( zeroorder + firstorder * (u - ul) + secondorder * (u - ul)^(2) )
end


"""

Theta, the anomaly for computing orbit averages
as a function of (a,e)

@IMPROVE, once happy with this version, remove ThetaRpRa

"""
function ThetaAE(ψ::Function,
                 dψ::Function,
                 d2ψ::Function,
                 d3ψ::Function,
                 u::Float64,
                 a::Float64,
                 e::Float64;
                 EDGE::Float64=0.01,
                 TOLECC::Float64=ELTOLECC,
                 f::Function=henon_f,
                 df::Function=henon_df,
                 d2f::Function=henon_d2f,
                 d3f::Function=henon_d3f,
                 d4f::Function=henon_d4f)

    # use the expanded approximation
    if ((1-abs(u))<EDGE)
        return ThetaExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e;TOLECC=TOLECC,f=f,d2f=d2f,d3f=d3f,d4f=d4f)

    # if not close to the boundary, can calculate as normal
    else
        E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)
        #r = ru(f,u,a,e)
        r = ru(u,a,e)

        # this can somehow be negative: do we need an extra check?
        denomsq = 2*(E - ψeff(ψ,r,L))

        if denomsq < 0.0
            # go back to the expansion -- or should we return 0.0?
            return ThetaExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e;TOLECC=TOLECC,f=f,d2f=d2f,d3f=d3f,d4f=d4f)
        end

        # do the standard return
        return a * e * df(u) / sqrt(denomsq)
    end

end


"""
numerical differentiation of Theta w.r.t. semimajor axis and eccentricity

"""
function ThetaAEdade(ψ::Function,
                     dψ::Function,
                     d2ψ::Function,
                     d3ψ::Function,
                     u::Float64,
                     a::Float64,
                     e::Float64;
                     EDGE::Float64=0.01,
                     TOLECC::Float64=ELTOLECC,
                     f::Function=henon_f,
                     df::Function=henon_df,
                     d2f::Function=henon_d2f,
                     d3f::Function=henon_d3f,
                     d4f::Function=henon_d4f,
                     da::Float64=1.0e-8,
                     de::Float64=1.0e-8)

    # derivative w.r.t. semimajor axis: always safe
    thHa = ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a+da,e,EDGE=EDGE,TOLECC=TOLECC,f=f,df=df,d2f=d2f,d3f=d3f,d4f=d4f)
    thLa = ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a   ,e,EDGE=EDGE,TOLECC=TOLECC,f=f,df=df,d2f=d2f,d3f=d3f,d4f=d4f)
    dThetada = (thHa-thLa)/da

    # derivative w.r.t. semimajor axis: safe unless too close to e=1
    # in that case, reverse de
    if e>(1-de)
        de *= -1
    end

    thHe = ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a,e+de,EDGE=EDGE,TOLECC=TOLECC,f=f,df=df,d2f=d2f,d3f=d3f,d4f=d4f)
    thLe = ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a,e   ,EDGE=EDGE,TOLECC=TOLECC,f=f,df=df,d2f=d2f,d3f=d3f,d4f=d4f)
    dThetade = (thHe-thLe)/de

    return dThetada,dThetade

end
