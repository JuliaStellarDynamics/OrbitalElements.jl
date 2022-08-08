"""
Ufunc.jl

collection of functions that are used to compute quantities when integrating over an orbit angle.

Specific to Henon mapping.

"""

"""ru(u,rp,ra)
mapping from u->r in Henon variables
"""
function ru(u::Float64,rp::Float64,ra::Float64)
    a,e = ae_from_rpra(rp,ra)
    fu = u*(3/2 - (u^2)/2)
    return a*(1+e*fu)
end

"""drdu(u,rp,ra)
first derivative of r(u)
"""
function drdu(u::Float64,rp::Float64,ra::Float64)
    a,e = ae_from_rpra(rp,ra)
    return (3/2)*a*e*(1-(u^2))
end

"""d2rdu(u,rp,ra)
second derivative of r(u)
"""
function d2rdu(u::Float64,rp::Float64,ra::Float64)
    a,e = ae_from_rpra(rp,ra)
    return -3*a*e*u
end

"""d3rdu(u,rp,ra)
third derivative of r(u)
"""
function d3rdu(u::Float64,rp::Float64,ra::Float64)
    a,e = ae_from_rpra(rp,ra)
    return -3*a*e
end

"""
all additional derivatives of r(u) are 0
"""

"""Q(ψ,dψ/dr,d²ψ/dr²,u,rp,ra)
Q \equiv 2(E-psi(r)) - L^2/r^2
uses mapping r(u)
"""
function Q(ψ::Function,
           dψ::Function,
           d2ψ::Function,
           u::Float64,
           rp::Float64,
           ra::Float64)

    E = E_from_rpra_pot(ψ,dψ,d2ψ,rp,ra)
    L = L_from_rpra_pot(ψ,dψ,d2ψ,rp,ra)

    r  = ru(u,rp,ra)

    if L==0
        return 2*(E-ψ(r))
    else
        return 2*(E-ψ(r)) - (L^2)/(r^2)
    end
end

"""expandTheta(ψ,dψ/dr,d²ψ/dr²,u,rp,ra[,VERBOSE])
expansion of Theta near u=pm1.
"""
function expandTheta(ψ::Function,
                     dψ::Function,
                     d2ψ::Function,
                     u::Float64,
                     rp::Float64,
                     ra::Float64;
                     VERBOSE::Int64=0,
                     SECONDORDER::Bool=true)

    # check for radial orbits
    L = L_from_rpra_pot(ψ,dψ,d2ψ,rp,ra)

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
    ddr = d2rdu(usign,rp,ra)
    dddr = d3rdu(usign,rp,ra)

    # if the sqrt is about to be zero...don't let it.
    #if ((ddr < 0) ⊻ (dQdrval < 0)) # XOR gate method, is this faster or slower?
    if (ddr*dQdrval < 0)

        if VERBOSE>0
            println("OrbitalElements/Ufunc.jl: Bad Theta expansion=",ddr,",",dQdrval," (u=",u,",rp=",rp,",ra=",ra,")")
        end

        u = -1.
        ddr = d2rdu(u,rp,ra)
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
function Theta(ψ::Function,
               dψ::Function,
               d2ψ::Function,
               u::Float64,
               rp::Float64,
               ra::Float64;
               EDGE::Float64=0.01,
               VERBOSE::Int64=0,
               SECONDORDER::Bool=true)

    L = L_from_rpra_pot(ψ,dψ,d2ψ,rp,ra)

    if (L==0)
        EDGE=1.e-10
    end

    if ((1-abs(u))<EDGE)
        # if near the edge, use the expanded approximations
        res = expandTheta(ψ,dψ,d2ψ,u,rp,ra)

    else # main branch, direct calculation

        Qval = Q(ψ,dψ,d2ψ,u,rp,ra)

        # check if Q is negative. If negative...revert back to the expanded values
        if (Qval<0)
            if VERBOSE>0
                println("OrbitalElements/Ufunc.jl: Bad Q=",Qval," (u=",u,",rp=",rp,",ra=",ra,")")
            end

            res = expandTheta(ψ,dψ,d2ψ,u,rp,ra)

        # if not negative, proceed along main branch
        else
            dr  = drdu(u,rp,ra)
            res = (1/sqrt(Qval))*dr
        end

    end

    return res
end

"""dQdr(ψ,dψ/dr,d²ψ/dr²,u,rp,ra)

Partial derivative of Q w.r.t. r

@IMPROVE, watch for switch when L==0

"""
function dQdr(ψ::Function,dψ::Function,d2ψ::Function,u::Float64,rp::Float64,ra::Float64)

    r   = ru(u,rp,ra)
    L = L_from_rpra_pot(ψ,dψ,d2ψ,rp,ra)

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

    r   = ru(u,rp,ra)
    L = L_from_rpra_pot(ψ,dψ,d2ψ,rp,ra)

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

    r   = ru(u,rp,ra)
    dr  = drdu(u,rp,ra)
    ddr = d2rdu(u,rp,ra)

    dQdu = dQdr(ψ,dψ,d2ψ,u,rp,ra)*drdu(u,rp,ra)

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
    r1   = ru(u,rp,ra)
    r2   = ru(u+eps,rp,ra)
    return (Qval2-Qval1)/(r2-r1)
end

"""dQdunumerical(ψ,dψ/dr,d²ψ/dr²,u,rp,ra[,eps])
numerical check version of dQ/du
"""
function dQdunumerical(ψ::Function,dψ::Function,d2ψ::Function,u::Float64,rp::Float64,ra::Float64,eps::Float64=0.001)
    Qval1 = Q(ψ,dψ,d2ψ,u,rp,ra)
    Qval2 = Q(ψ,dψ,d2ψ,u+eps,rp,ra)
    return (Qval2-Qval1)/(eps)
end

"""dThetadunumerical(ψ,dψ/dr,d²ψ/dr²,u,rp,ra)
numerical check version of dTheta/du
"""
function dThetadunumerical(ψ::Function,dψ::Function,d2ψ::Function,u::Float64,rp::Float64,ra::Float64)
    eps = 0.001
    Tval1 = Theta(ψ,dψ,d2ψ,u,rp,ra,0.)
    Tval2 = Theta(ψ,dψ,d2ψ,u+eps,rp,ra,0.)
    return (Tval2-Tval1)/(eps)
end
