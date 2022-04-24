"""
Ufunc.jl

collection of functions that are used to compute quantities when integrating over an orbit angle.

Specific to Henon mapping.

"""

"""ru

mapping from u->r in Henon variables

"""
function ru(u::Float64,rperi::Float64,rapo::Float64)
    a,e = (rperi+rapo)/2,(rapo-rperi)/(rapo+rperi)
    fu = u*(3/2 - (u^2)/2)
    return a*(1+e*fu)
end

"""drdu

first derivative of r(u)

"""
function drdu(u::Float64,rperi::Float64,rapo::Float64)
    a,e = (rperi+rapo)/2,(rapo-rperi)/(rapo+rperi)

    return (3/2)*a*e*(1-(u^2))
end

"""ddrddu

second derivative of r(u)

"""
function ddrddu(u::Float64,rperi::Float64,rapo::Float64)
    a,e = (rperi+rapo)/2,(rapo-rperi)/(rapo+rperi)

    return -3*a*e*u
end

"""Q(u,rp,ra)

Q \equiv 2(E-psi(r)) - L^2/r^2

uses mapping r(u)

"""
function Q(potential::Function,dpotential::Function,ddpotential::Function,u::Float64,rp::Float64,ra::Float64)

    E = E_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
    L = L_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)

    r  = ru(u,rp,ra)
    #println(E," ",L," ",r," ",potential(r))

    if L==0
        return 2*(E-potential(r))# - (L^2)/(r^2)
    else
        return 2*(E-potential(r)) - (L^2)/(r^2)
    end
end

"""expandThetaNeg

expansion of Theta near -1.

"""
function expandThetaNeg(potential::Function,dpotential::Function,ddpotential::Function,u::Float64,rp::Float64,ra::Float64)

    # check for radial orbits
    L = L_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)

    if (L==0)
        return 0.
    end

    ddr = ddrddu(u,rp,ra)
    dQdrval = dQdr(potential,dpotential,ddpotential,u,rp,ra)

    # if the sqrt is about to be zero...don't let it.
    if ((ddr < 0) | (dQdrval < 0))
        println("OrbitalElements/Ufunc.jl: Bad T-expansion=",ddr,",",dQdrval," (u=",u,",rp=",rp,",ra=",ra,")")
        u = -1.
        ddr = ddrddu(u,rp,ra)
        dQdrval = dQdr(potential,dpotential,ddpotential,u,rp,ra)
    end

    return sqrt(2*ddr/dQdrval)
end

"""expandThetaPos

expansion of Theta near +1.

"""
function expandThetaPos(potential::Function,dpotential::Function,ddpotential::Function,u::Float64,rp::Float64,ra::Float64)
    ddr = ddrddu(u,rp,ra)
    dQdrval = dQdr(potential,dpotential,ddpotential,u,rp,ra)

    # if the sqrt is about to be zero...don't let it.
    if ((ddr < 0) | (dQdrval < 0))
        println("OrbitalElements/Ufunc.jl: Bad T-expansion=",ddr,",",dQdrval," (u=",u,",rp=",rp,",ra=",ra,")")
        u = 1.
        ddr = ddrddu(u,rp,ra)
        dQdrval = dQdr(potential,dpotential,ddpotential,u,rp,ra)
    end

    return sqrt(abs(2*ddr/dQdrval))
end


"""Theta

``Theta(s) \\equiv (1/vr)*(dr/ds) where s\\in[-1,1]``

ulim sets the switch to the expanded approximation near the boundaries
for radial orbits, the pericentre limit is overridden and set to 0.
"""
function Theta(potential::Function,dpotential::Function,ddpotential::Function,u::Float64,rp::Float64,ra::Float64,ulim::Float64=0.01)

    L = L_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
    if (L==0)
        ulim=1.e-10
    end

    if ((1-abs(u))<ulim)# use the expanded approximations

        if u>0
            res = expandThetaPos(potential,dpotential,ddpotential,u,rp,ra)
        else
            res = expandThetaNeg(potential,dpotential,ddpotential,u,rp,ra)
        end

    else # main branch, direct calculation

        Qval = Q(potential,dpotential,ddpotential,u,rp,ra)

        # check if Q is negative. If negative...revert back to the expanded values
        if (Qval<0)
            println("OrbitalElements/Ufunc.jl: Bad Q=",Qval," (u=",u,",rp=",rp,",ra=",ra,")")

            if (u>0)
                res = expandThetaPos(potential,dpotential,ddpotential,u,rp,ra)
            else
                res = expandThetaNeg(potential,dpotential,ddpotential,u,rp,ra)
            end

        # if not negative, proceed along main branch
        else
            dr = drdu(u,rp,ra)
            res = (1/sqrt(Qval))*dr
        end


    end

    return res
end

"""dQdr

Partial derivative of Q w.r.t. r

@IMPROVE, watch for switch when L==0

"""
function dQdr(potential::Function,dpotential::Function,ddpotential::Function,u::Float64,rp::Float64,ra::Float64)

    r   = ru(u,rp,ra)
    L = L_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)

    if (L==0)
        return max(0.,-2*dpotential(r)-2*ddpotential(r)*r)
    else
        return  ((2*(L^2)/(r^3)) - 2*dpotential(r))
    end

end

"""dThetadu

partial derivative of Theta w.r.t. u.

"""
function dThetadu(potential::Function,dpotential::Function,ddpotential::Function,u::Float64,rp::Float64,ra::Float64)

    Qval = Q(potential,dpotential,ddpotential,u,rp,ra)
    #E = E_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
    #L = L_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)

    r   = ru(u,rp,ra)
    dr  = drdu(u,rp,ra)
    ddr = ddrddu(u,rp,ra)

    dQdu = dQdr(potential,dpotential,ddpotential,u,rp,ra)*drdu(u,rp,ra)

    return (Qval^(-1/2))*ddr + dr*((-1/2)*dQdu*(Qval^(-3/2)))

end

"""dQdrnumerical

numerical check version of dQ/dr

"""
function dQdrnumerical(potential::Function,dpotential::Function,ddpotential::Function,u::Float64,rp::Float64,ra::Float64,eps::Float64=0.001)
    Qval1 = Q(potential,dpotential,ddpotential,u,rp,ra)
    Qval2 = Q(potential,dpotential,ddpotential,u+eps,rp,ra)
    r1   = ru(u,rp,ra)
    r2   = ru(u+eps,rp,ra)
    return (Qval2-Qval1)/(r2-r1)
end

"""dQdunumerical

numerical check version of dQ/du

"""
function dQdunumerical(potential::Function,dpotential::Function,ddpotential::Function,u::Float64,rp::Float64,ra::Float64,eps::Float64=0.001)
    Qval1 = Q(potential,dpotential,ddpotential,u,rp,ra)
    Qval2 = Q(potential,dpotential,ddpotential,u+eps,rp,ra)
    return (Qval2-Qval1)/(eps)
end

"""dThetadunumerical

numerical check version of dTheta/du

"""
function dThetadunumerical(potential::Function,dpotential::Function,ddpotential::Function,u::Float64,rp::Float64,ra::Float64)
    eps = 0.001
    Tval1 = Theta(potential,dpotential,ddpotential,u,rp,ra,0.)
    Tval2 = Theta(potential,dpotential,ddpotential,u+eps,rp,ra,0.)
    return (Tval2-Tval1)/(eps)
end
