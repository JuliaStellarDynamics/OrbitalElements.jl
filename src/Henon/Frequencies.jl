

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
                                   r_apo::Float64,
                                   r_peri::Float64,
                                   E::Float64,
                                   L::Float64,
                                   action::Bool=false,
                                   NINT::Int64=64)



    # set the integration width
    integration_distance = 2
    du = integration_distance/NINT;

    # define some auxilliaries
    # @IMPROVE: offload these to OrbitDefinitions.jl
    ap  = 0.5*(r_apo + r_peri); # a
    am  = 0.5*(r_apo - r_peri); # ae
    ecc = (r_apo-r_peri)/(r_apo+r_peri)#am/ap
    sp  = ap/(r_apo*r_peri);
    sm  = am/(r_apo*r_peri);

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

        dr = (3/4)*(r_apo-r_peri)*(1-(u^2))

        ur = potential(r)

        tmp = 2(E-ur) - (L*L)/(r*r)

        s    = r/(r_apo*r_peri);
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
        r_peri::Float64,r_apo::Float64,action::Bool=false,TOLECC::Float64=0.01,verbose::Int64=0)

    E = E_from_rpra_pot(potential,dpotential,ddpotential,r_peri,r_apo)
    J = L_from_rpra_pot(potential,dpotential,ddpotential,r_peri,r_apo)


    # check the tolerance
    a,ecc = ae_from_rpra(r_peri,r_apo)

    # don't go into the loop if circular
    if ecc<TOLECC
        return Omega1_circular(dpotential,ddpotential,a),Omega2_circular(dpotential,a)
    end

    # don't go into the loop if radial
    if (1-ecc)<TOLECC
        #print("Too radial!\n")
        freq1,freq2 = henon_anomaly_frequencies(potential,r_apo,1.e-10,E,J)
        return freq1,freq2
    end

    # go to the frequency calculation
    if action
        freq1,freq2,action1 = henon_anomaly_frequencies(potential,r_apo,r_peri,E,J,true)
        return freq1,freq2,action1
    else
        freq1,freq2 = henon_anomaly_frequencies(potential,r_apo,r_peri,E,J,false)
        return freq1,freq2
    end
end

"""compute_frequences_henon_ae(ψ,dψ/dr,d²ψ/dr²,rp,ra[,TOLECC,verbose])
"""
function compute_frequencies_henon_ae(potential::Function,dpotential::Function,ddpotential::Function,
                                      a::Float64,ecc::Float64,action::Bool=false,TOLECC::Float64=0.01,verbose::Int64=0)
    #=compute_frequencies_henon_ae

    wrapper to compute the frequencies for an orbit specified by (a,e)

    =#

    # @IMPROVE: go exactly to 1.0 in eccentricity
    # if too radial, don't let J go to zero
    if (1-ecc)<0.01*TOLECC
        r_peri,r_apo = rpra_from_ae(a,0.999999)
    else
        r_peri,r_apo = rpra_from_ae(a,ecc)
    end


    E = E_from_rpra_pot(potential,dpotential,ddpotential,r_peri,r_apo)
    J = L_from_rpra_pot(potential,dpotential,ddpotential,r_peri,r_apo)

    if verbose>0
        print("E/J ",E," ",J,"\n")
    end

    # don't go into the loop if circular
    if ecc<TOLECC
        return Omega1_circular(dpotential,ddpotential,a),Omega2_circular(dpotential,a)
    end

    # if radial, we should exploit freq1 = 2*freq2 somehow.

    # go to the frequency calculation
    if action
        freq1,freq2,action1 = henon_anomaly_frequencies(potential,r_apo,r_peri,E,J,true)
        return freq1,freq2,action1
    else
        freq1,freq2 = henon_anomaly_frequencies(potential,r_apo,r_peri,E,J,false)
        return freq1,freq2
    end

end
