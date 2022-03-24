

function u_potential_rpra(potential::Function,dpotential::Function,ddpotential::Function,
                          u::Float64,r_peri::Float64,r_apo::Float64,EDGETOL::Float64=1.e-5)
    #= Q(s) \equiv (1/vr)*(dr/ds) where s\in[-1,1].


    =#

    J = L_from_rpra_pot(potential,dpotential,ddpotential,r_peri,r_apo,0.5*EDGETOL)

    edge = 0.

    if u<-1+EDGETOL
        edge = g_minus1(potential,dpotential,ddpotential,a,ecc)
        return edge
    end

    if u>1-EDGETOL
        edge = g_plus1(potential,dpotential,ddpotential,a,ecc)
        return edge
    end

    fu = u*(3/2 - u*u/2)
    r  = a*(1+ecc*fu)

    # clearly when dr = +/- 1 this will fail! dr -> 0.
    # maybe do an expansion?
    #dr = (3/4)*(2a*ecc)*(1-(u^2))
    dr = drdu(u,a,ecc)

    #print(r," ",dr," ",fu,"\n")

    ur = potential(r)
    #print("UR=",ur,"\n")

    E = E_from_rpra_pot(potential,dpotential,ddpotential,r_peri,r_apo,0.5*EDGETOL)

    # safe when sufficiently far from u \pm 1
    tmp = 2(E-ur) - (J*J)/(r*r)


    best = dr/sqrt(maximum([1.e-16,tmp]))

    if isnan(best)
        return 0.0
    else
        return best
    end

end



function u_potential_aecc(potential::Function,dpotential::Function,ddpotential::Function,
                          u::Float64,a::Float64,ecc::Float64,EDGETOL::Float64=1.e-5)
    #= Q(s) \equiv (1/vr)*(dr/ds) where s\in[-1,1].


    =#

    r_peri,r_apo = rpra_from_ae(a,ecc)
    J = L_from_rpra_pot(potential,dpotential,ddpotential,r_peri,r_apo,0.5*EDGETOL)

    edge = 0.

    if u<-1+EDGETOL
        edge = g_minus1(potential,dpotential,ddpotential,a,ecc)
        return edge
    end

    if u>1-EDGETOL
        edge = g_plus1(potential,dpotential,ddpotential,a,ecc)
        return edge
    end

    fu = u*(3/2 - u*u/2)
    r  = a*(1+ecc*fu)

    # clearly when dr = +/- 1 this will fail! dr -> 0.
    # maybe do an expansion?
    #dr = (3/4)*(2a*ecc)*(1-(u^2))
    dr = drdu(u,a,ecc)

    #print(r," ",dr," ",fu,"\n")

    ur = potential(r)
    #print("UR=",ur,"\n")

    E = E_from_rpra_pot(potential,dpotential,ddpotential,r_peri,r_apo,0.5*EDGETOL)

    # safe when sufficiently far from u \pm 1
    tmp = 2(E-ur) - (J*J)/(r*r)


    best = dr/sqrt(maximum([1.e-16,tmp]))

    # do an interpolation:
    #if (edge>0) & (u>0)
    #    dr1 = 1 - u
    #    dr2 = u - 0.8
    #    #print(dr1," ",dr2)
    #    return best*dr2 + edge*dr1
    #end

    if isnan(best)
        return 0.0
    else
        return best
    end

end


function drdu(u::Float64,a::Float64,ecc::Float64,TOL::Float64=0.01)
    #=(dr/ds), the Jacobian of the change of variables from r to s.

    can include expanded values so (dr/ds) goes to zero at a faster rate than (1/vr)

    =#
    dr = (3/4)*(2a*ecc)*(1-(u^2))
    ddr = -3a*ecc*u

    if (1 - u^2)<TOL
        # expand from value
        if u > 0 # u near 1
            # u + h, h is negative, u is positive
            h = u-1
            dr = h*ddr
        else # u near -1
            # u + h, h is positive, u is negative
            h = u+1
            dr = h*ddr
        end
    end

    return dr

end



function dQdr(dpotential::Function,r::Float64,L::Float64)
    #=the derivative of Q = (1/vr)*(dr/ds) w.r.t. r

    =#

    return (2(L^2)/(r^3)) - 2*dpotential(r)
end

function ddQddr(ddpotential::Function,r::Float64,L::Float64)
    #=the second derivative of Q = (1/vr)*(dr/ds) w.r.t. r

    =#
    return (-6(L^2)/(r^4)) - 2*ddpotential(r)
end



function g_minus1(potential::Function,dpotential::Function,ddpotential::Function,
                 a::Float64,ecc::Float64)
    #= The edge expansion from the change of variables.
    See Henon 1971, eq. 52


    =#

    if ecc<0.00001
        ecc = 0.00001
    end

    if ecc>0.99999999
        return 0.
    end

    rperi,rapo = rpra_from_ae(a,ecc)
    L = L_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo,0.0000005)

    dq = dQdr(dpotential,rperi,L)
    #print(dq,"\n")

    return sqrt(3(rapo-rperi)/dq)
end

function g_plus1(potential::Function,dpotential::Function,ddpotential::Function,
                 a::Float64,ecc::Float64)

    if ecc<0.00001
        ecc = 0.00001
    end

    if ecc>0.9999999
        ecc = 0.9999999
    end

    rperi,rapo = rpra_from_ae(a,ecc)
    L = L_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo,0.0000005)

    dq = dQdr(dpotential,rapo,L)
    #print(dq,"\n")

    return sqrt(maximum([-3(rapo-rperi)/dq,0]))
end
