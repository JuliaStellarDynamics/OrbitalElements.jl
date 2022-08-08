

function henon_theta_frequencies(potential::Function,
                                 dpotential::Function,
                                 ddpotential::Function,
                                 rperi::Float64,
                                 rapo::Float64;
                                 action::Bool=false,
                                 NINT::Int64=1000,
                                 EDGE::Float64=0.02)

    # using theta calculations to compute frequencies: leans heavily on Theta from Ufunc.jl
    # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

    # currently using Simpson's 1/3 rule for integration: requires that NINT be even.
    # @IMPROVE: we could use a better integration scheme

    # set integrations to zero: omega1, omega2, action1
    accum1,accum2,accum3 = 0.0,0.0,0.0

    # build composite Simpson's rule sample points
    x = LinRange(-1.0,1.0,NINT)

    # start with the endpoint weights
    accum1 += Theta(potential,dpotential,ddpotential,x[1],rperi,rapo,EDGE)
    accum1 += Theta(potential,dpotential,ddpotential,x[NINT],rperi,rapo,EDGE)

    accum2 += Theta(potential,dpotential,ddpotential,x[1],rperi,rapo,EDGE)/(ru(x[NINT],rperi,rapo)^2)
    accum2 += Theta(potential,dpotential,ddpotential,x[NINT],rperi,rapo,EDGE)/(ru(x[NINT],rperi,rapo)^2)

    accum3 += Theta(potential,dpotential,ddpotential,x[1],rperi,rapo,EDGE)*Q(potential,dpotential,ddpotential,x[1],rperi,rapo)
    accum3 += Theta(potential,dpotential,ddpotential,x[NINT],rperi,rapo,EDGE)*Q(potential,dpotential,ddpotential,x[NINT],rperi,rapo)


    # now accumulate first set of interior points
    for j=1:(NINT÷2)
        accum1 += 4 * Theta(potential,dpotential,ddpotential,x[2j-1],rperi,rapo,EDGE)
        accum2 += 4 * Theta(potential,dpotential,ddpotential,x[2j-1],rperi,rapo,EDGE)/(ru(x[2j-1],rperi,rapo)^2)
        accum3 += 4 * Theta(potential,dpotential,ddpotential,x[2j-1],rperi,rapo,EDGE) * Q(potential,dpotential,ddpotential,x[2j-1],rperi,rapo)
    end

    # now accumulate second set of interior points
    for j=1:((NINT÷2) - 1)
        accum1 += 2 * Theta(potential,dpotential,ddpotential,x[2j],rperi,rapo,EDGE)
        accum2 += 2 * Theta(potential,dpotential,ddpotential,x[2j],rperi,rapo,EDGE)/(ru(x[2j],rperi,rapo)^2)
        accum3 += 2 * Theta(potential,dpotential,ddpotential,x[2j],rperi,rapo,EDGE) * Q(potential,dpotential,ddpotential,x[2j],rperi,rapo)
    end

    # complete the integration with weighting
    h = 2/(3NINT)
    accum1 *= h
    accum2 *= h
    accum3 *= h

    #return the values
    omega1inv = (1/pi)*accum1
    omega1    = 1/omega1inv
    omega2    = L_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo)*accum2*(1/pi)*omega1
    actionj   = accum3/pi

    if action
        return omega1,omega2,actionj
    else
        return omega1,omega2
    end

end




"""the henon anomaly increment
"""
function henon_f(u::Float64)
    u*(1.5 - 0.5*u^(2))
end

"""the derivative of the henon anomaly increment
"""
function henon_dfdu(u::Float64)
    1.5*(1.0 - u^(2))
end



"""henon_anomaly_frequencies

use the henon anomaly mapping to compute orbit frequencies

@IMPROVE, handle negative sqrts gracefully
@IMPROVE, add switch for action calculation
@IMPROVE, offload transformations to OrbitDefinitions.jl
@IMPROVE, disallow any frequency overshoots with known boundaries

"""
function henon_anomaly_frequencies(potential::Function,
                                   rapo::Float64,
                                   rperi::Float64,
                                   E::Float64,
                                   L::Float64;
                                   action::Bool=false,
                                   NINT::Int64=64)



    # set the integration width
    integration_distance = 2
    du = integration_distance/NINT;

    # define some auxilliaries
    # @IMPROVE: offload these to OrbitDefinitions.jl
    ap  = 0.5*(rapo + rperi); # a
    am  = 0.5*(rapo - rperi); # ae
    ecc = (rapo-rperi)/(rapo+rperi)#am/ap
    sp  = ap/(rapo*rperi);
    sm  = am/(rapo*rperi);

    # set the accumulators to zero
    accum0 = 0.0
    accum1 = 0.0
    accum2 = 0.0

    # keep track of how many steps we actually use
    uNINT = 0.0

    # proceed as centred rectangle integration: starting point
    u = 0.5*(du-integration_distance)

    for i=1:NINT

        fu = u*(3/2 - u*u/2)
        r  = ap*(1+ecc*fu)

        dr = (3/4)*(rapo-rperi)*(1-(u^2))

        ur = potential(r)

        tmp = 2(E-ur) - (L*L)/(r*r)

        s    = r/(rapo*rperi);
        ur1  = potential(1.0/s)
        tmp2 = 2*(E-ur1) - (L*L*s*s)



        # @IMPROVE: this is a block for invalid accumulations, but this check could be smarter
        if (tmp>0) & (tmp2>0)
            accum0 += dr * sqrt(tmp);
            accum1 += dr / sqrt(tmp);
            accum2 += dr/sqrt(tmp2);
            uNINT += 1.0
        end

        # advance the counter
        u += du

    end


    #freq1 = integration_distance/(accum1*dt);
    freq1 = (pi/2)*uNINT/(accum1)
    freq2 = freq1/(pi/2) * L * (sm/am) * accum2 / uNINT;

    # @IMPROVE
    # note that we could also return the actions if we wanted:
    action1 = accum0/((pi/2)*uNINT)
    #action2 = jj;

    # @IMPROVE: we may want to force never allowing an overshoot (i.e. freq1 is capped at 1 no matter what)
    if action
        return freq1,freq2,action1
    else
        return freq1,freq2
    end

end


"""compute_frequences_henon(ψ,dψ/dr,d²ψ/dr²,rp,ra[,TOLECC,verbose])
"""
function compute_frequencies_henon(potential::Function,dpotential::Function,ddpotential::Function,
        rperi::Float64,rapo::Float64;action::Bool=false,TOLECC::Float64=0.01,verbose::Int64=0,NINT::Int64=32)

    E = E_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo)
    J = L_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo)


    # check the tolerance
    a,ecc = ae_from_rpra(rperi,rapo)

    # don't go into the loop if circular
    if ecc<TOLECC
        return Omega1_circular(dpotential,ddpotential,a),Omega2_circular(dpotential,a)
    end

    # don't go into the loop if radial
    if (1-ecc)<TOLECC
        rperi = 1.e-10
        if action
            E = E_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo)
            J = L_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo)
            freq1,freq2,action1 = henon_anomaly_frequencies(potential,rapo,rperi,E,J,action=true,NINT=NINT)
            return freq1,freq2
        else
            E = E_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo)
            J = L_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo)
            freq1,freq2 = henon_anomaly_frequencies(potential,rapo,rperi,E,J,action=false,NINT=NINT)
            return freq1,freq2
        end
    end

    # go to the frequency calculation
    if action
        freq1,freq2,action1 = henon_anomaly_frequencies(potential,rapo,rperi,E,J,action=true,NINT=NINT)
        return freq1,freq2,action1
    else
        freq1,freq2 = henon_anomaly_frequencies(potential,rapo,rperi,E,J,action=false,NINT=NINT)
        return freq1,freq2
    end
end

"""compute_frequences_henon_ae(ψ,dψ/dr,d²ψ/dr²,rp,ra[,TOLECC,verbose])
"""
function compute_frequencies_henon_ae(potential::Function,dpotential::Function,ddpotential::Function,
                                      a::Float64,ecc::Float64;action::Bool=false,TOLECC::Float64=0.01,verbose::Int64=0,NINT=32)
    #=compute_frequencies_henon_ae

    wrapper to compute the frequencies for an orbit specified by (a,e)

    =#

    # @IMPROVE: go exactly to 1.0 in eccentricity
    # if too radial, don't let J go to zero
    if (1-ecc)<0.01*TOLECC
        rperi,rapo = rpra_from_ae(a,0.999999)
    else
        rperi,rapo = rpra_from_ae(a,ecc)
    end


    E = E_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo)
    J = L_from_rpra_pot(potential,dpotential,ddpotential,rperi,rapo)

    if verbose>2
        print("E/J ",E," ",J,"\n")
    end

    # don't go into the loop if too circular: use an approximation
    if ecc<TOLECC
        if action
            return Omega1_circular(dpotential,ddpotential,a),Omega2_circular(dpotential,a),0.0
        else
            return Omega1_circular(dpotential,ddpotential,a),Omega2_circular(dpotential,a)
        end
    end

    # if radial, we should exploit freq1 = 2*freq2 somehow.

    # go to the frequency calculation
    if action
        freq1,freq2,action1 = henon_anomaly_frequencies(potential,rapo,rperi,E,J,action=true,NINT=NINT)
        return freq1,freq2,action1
    else
        freq1,freq2 = henon_anomaly_frequencies(potential,rapo,rperi,E,J,action=false,NINT=NINT)
        return freq1,freq2
    end

end
