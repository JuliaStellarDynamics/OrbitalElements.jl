




"""HenonThetaFrequencies(ψ,dψ,d2ψ,rp,ra[,action=false,NINT=32,EDGE=0.03,TOLECC=0.001])

use the defined function Theta(u) to compute frequency integrals
"""
function HenonThetaFrequencies(ψ::Function,
                               dψ::Function,
                               d2ψ::Function,
                               rperi::Float64,
                               rapo::Float64;
                               action::Bool=false,
                               NINT::Int64=32,
                               EDGE::Float64=0.03,
                               TOLECC::Float64=0.001)


    # check the eccentricity tolerances
    a,e = ae_from_rpra(rperi,rapo)

    if e<TOLECC
        # drop into circular frequency calculation
        if action
            return Omega1_circular(dψ,d2ψ,a),Omega2_circular(dψ,a),0.0
        else
            return Omega1_circular(dψ,d2ψ,a),Omega2_circular(dψ,a)
        end

    else
        # make a hard barrier for orbit calculation to avoid too radial orbits
        if e>(1-TOLECC)
            rperi,rapo = rpra_from_ae(a,1-TOLECC)
        end

        # standard case, compute the hard integrals


        # using theta calculations to compute frequencies: leans heavily on Theta from Ufunc.jl
        # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

        # currently using Simpson's 1/3 rule for integration: requires that NINT be even.
        # @IMPROVE: we could use a better integration scheme

        function u3func(u::Float64)
            # push integration forward on three different quantities: Theta(u),Theta(u)/r^2(u),Theta(u)*vr(u)

            th = Theta(ψ,dψ,d2ψ,u,rperi,rapo,EDGE=EDGE)

            return (th,
                    th/(ru(u,rperi,rapo)^2),
                    th*(Q(ψ,dψ,d2ψ,u,rperi,rapo)))

        end

        accum1,accum2,accum3 = UnitarySimpsonIntegration(u3func,K_O=NINT)

        #return the values
        omega1inv = (1/pi)*accum1
        omega1    = 1/omega1inv
        actionj   = (1/pi)*accum3

        # be careful with omega2 if near radial: use analytic relation
        if e>(1-TOLECC)
            omega2 = 0.5*omega1
        else
            omega2 = LFromRpRa(ψ,dψ,d2ψ,rperi,rapo)*accum2*(1/pi)*omega1
        end

        if action
            return omega1,omega2,actionj
        else
            return omega1,omega2
        end

    end # switches for orbits who are too radial or circular

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
function henon_anomaly_frequencies(ψ::Function,
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

        ur = ψ(r)

        tmp = 2(E-ur) - (L*L)/(r*r)

        s    = r/(rapo*rperi);
        ur1  = ψ(1.0/s)
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
function compute_frequencies_henon(ψ::Function,dψ::Function,d2ψ::Function,
        rperi::Float64,rapo::Float64;action::Bool=false,TOLECC::Float64=0.01,verbose::Int64=0,NINT::Int64=32)

    E = EFromRpRa(ψ,dψ,d2ψ,rperi,rapo)
    J = LFromRpRa(ψ,dψ,d2ψ,rperi,rapo)


    # check the tolerance
    a,ecc = ae_from_rpra(rperi,rapo)

    # don't go into the loop if circular
    if ecc<TOLECC
        return Omega1_circular(dψ,d2ψ,a),Omega2_circular(dψ,a)
    end

    # don't go into the loop if radial
    if (1-ecc)<TOLECC
        rperi = 1.e-10
        if action
            E = EFromRpRa(ψ,dψ,d2ψ,rperi,rapo)
            J = LFromRpRa(ψ,dψ,d2ψ,rperi,rapo)
            freq1,freq2,action1 = henon_anomaly_frequencies(ψ,rapo,rperi,E,J,action=true,NINT=NINT)
            return freq1,freq2
        else
            E = EFromRpRa(ψ,dψ,d2ψ,rperi,rapo)
            J = LFromRpRa(ψ,dψ,d2ψ,rperi,rapo)
            freq1,freq2 = henon_anomaly_frequencies(ψ,rapo,rperi,E,J,action=false,NINT=NINT)
            return freq1,freq2
        end
    end

    # go to the frequency calculation
    if action
        freq1,freq2,action1 = henon_anomaly_frequencies(ψ,rapo,rperi,E,J,action=true,NINT=NINT)
        return freq1,freq2,action1
    else
        freq1,freq2 = henon_anomaly_frequencies(ψ,rapo,rperi,E,J,action=false,NINT=NINT)
        return freq1,freq2
    end
end

"""compute_frequences_henon_ae(ψ,dψ/dr,d²ψ/dr²,rp,ra[,TOLECC,verbose])
"""
function compute_frequencies_henon_ae(ψ::Function,dψ::Function,d2ψ::Function,
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


    E = EFromRpRa(ψ,dψ,d2ψ,rperi,rapo)
    J = LFromRpRa(ψ,dψ,d2ψ,rperi,rapo)

    if verbose>2
        print("E/J ",E," ",J,"\n")
    end

    # don't go into the loop if too circular: use an approximation
    if ecc<TOLECC
        if action
            return Omega1_circular(dψ,d2ψ,a),Omega2_circular(dψ,a),0.0
        else
            return Omega1_circular(dψ,d2ψ,a),Omega2_circular(dψ,a)
        end
    end

    # if radial, we should exploit freq1 = 2*freq2 somehow.

    # go to the frequency calculation
    if action
        freq1,freq2,action1 = henon_anomaly_frequencies(ψ,rapo,rperi,E,J,action=true,NINT=NINT)
        return freq1,freq2,action1
    else
        freq1,freq2 = henon_anomaly_frequencies(ψ,rapo,rperi,E,J,action=false,NINT=NINT)
        return freq1,freq2
    end

end
