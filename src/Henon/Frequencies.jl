
function henon_anomaly_frequencies(potential::Function,r_apo::Float64,r_peri::Float64,ee::Float64,jj::Float64)
    # set the integration width
    FRECS = 64
    integration_distance = 2
    du = integration_distance/FRECS;

    # define some auxilliaries
    ap  = 0.5*(r_apo + r_peri); # a
    am  = 0.5*(r_apo - r_peri); # ae
    ecc = (r_apo-r_peri)/(r_apo+r_peri)#am/ap
    sp  = ap/(r_apo*r_peri);
    sm  = am/(r_apo*r_peri);

    # set the accumulators to zero
    accum0 = 0.0;
    accum1 = 0.0;
    accum2 = 0.0;
    uFREQS = 0.0

    # proceed as centred rectangle integration: starting point
    u = 0.5*(du-integration_distance)

    for i=1:FRECS

        #print(i,' ',u,'\n')

        fu = u*(3/2 - u*u/2)
        r  = ap*(1+ecc*fu)
        #print(r," ")
        dr = (3/4)*(r_apo-r_peri)*(1-(u^2))
        #print(dr," ")

        ur = potential(r)
        #print(ur," ")

        tmp = 2(ee-ur) - (jj*jj)/(r*r)

        s    = r/(r_apo*r_peri);
        ur1  = potential(1.0/s)
        tmp2 = 2*(ee-ur1) - (jj*jj*s*s)
        #print(tmp," ")
        accum0 += dr * tmp;
        if (tmp>0) & (tmp2>0)
            accum1 += dr / sqrt(tmp);
            uFREQS += 1.0


            #if tmp>0
            accum2 += dr/sqrt(tmp2);

        end
        #print(accum1,"\n")


        #end

        # advance the counter
        u += du

    end


    #freq1 = integration_distance/(accum1*dt);
    freq1 = (pi/2)*uFREQS/(accum1)
    freq2 = freq1/(pi/2) * jj * (sm/am) * accum2 / uFREQS;

    # note that we could also compute the actions if we wanted:
    #action1 = am*accum0*dt/pi;
    #action2 = jj;

    # we may want to force never allowing an overshoot (i.e. freq1 is capped at 1 no matter what)
    return freq1,freq2

end



function compute_frequencies_henon(potential::Function,dpotential::Function,ddpotential::Function,
        r_peri::Float64,r_apo::Float64,
    TOLECC::Float64=0.01)

    E = E_from_rpra_pot(potential,dpotential,ddpotential,r_peri,r_apo)
    J = L_from_rpra_pot(potential,dpotential,ddpotential,r_peri,r_apo)

    #print("E/J ",E," ",J,"\n")


    # check the tolerance
    a,ecc = ae_from_rpra(r_peri,r_apo)

    # don't go into the loop if circular
    if ecc<TOLECC
        return Omega1_circular(dpotential,ddpotential,a),Omega2_circular(dpotential,a)
    end

    # don't go into the loop if radial
    if (1-ecc)<TOLECC
        print("Too radial!\n")
        freq1,freq2 = henon_anomaly_frequencies(potential,r_apo,1.e-10,E,J)
        return freq1,freq2
    end


    # go to the frequency calculation
    freq1,freq2 = henon_anomaly_frequencies(potential,r_apo,r_peri,E,J)

    return freq1,freq2
end

function compute_frequencies_henon_ae(potential::Function,dpotential::Function,ddpotential::Function,
                                      a::Float64,ecc::Float64,TOLECC::Float64=0.01,verbose::Int64=0)

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

    # go to the safe frequency calculation integration
    freq1,freq2 = henon_anomaly_frequencies(potential,r_apo,r_peri,E,J)
    return freq1,freq2

end
