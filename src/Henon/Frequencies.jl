

"""HenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e,[,action=false,NINT=32,EDGE=0.03,TOLECC=0.001])

use the defined function Theta(u) to compute frequency integrals
"""
function HenonThetaFrequenciesAE(ψ::Function,
                                 dψ::Function,
                                 d2ψ::Function,
                                 d3ψ::Function,
                                 a::Float64,
                                 e::Float64;
                                 action::Bool=false,
                                 NINT::Int64=32,
                                 EDGE::Float64=0.03,
                                 TOLECC::Float64=0.001)




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

        rperi,rapo = rpra_from_ae(a,e)

        function u3func(u::Float64)
            # push integration forward on three different quantities: Theta(u),Theta(u)/r^2(u),Theta(u)*vr(u)

            th = ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)

            return (th,
                    th/(ruAE(u,a,e)^2),
                    th*(Q(ψ,dψ,d2ψ,u,rperi,rapo)))

        end

        accum = UnitarySimpsonIntegration(u3func,K_O=NINT)

        #return the values
        Ω1inv = (1/pi)*accum[1]
        Ω1    = 1/Ω1inv
        actionj   = (1/pi)*accum[3]

        # be careful with Ω2 if near radial: use analytic relation
        if e>(1-TOLECC)
            Ω2 = 0.5*Ω1
        else
            Ω2 = LFromAE(ψ,dψ,d2ψ,d3ψ,a,e)*accum[2]*(1/pi)*Ω1
        end

        if action
            return Ω1,Ω2,actionj
        else
            return Ω1,Ω2
        end

    end # switches for orbits who are too radial or circular

end


"""DHenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,[da,de,action=false,NINT=32,EDGE=0.03,TOLECC=0.001])

use the defined function Theta(u) to compute frequency integrals
AND DERIVATIVES

"""
function DHenonThetaFrequenciesAE(ψ::Function,
                                  dψ::Function,
                                  d2ψ::Function,
                                  d3ψ::Function,
                                  d4ψ::Function,
                                  a::Float64,
                                  e::Float64;
                                  da::Float64=1.0e-6,
                                  de::Float64=1.0e-6,
                                  NINT::Int64=32,
                                  EDGE::Float64=0.01,
                                  TOLECC::Float64=0.001,
                                  Omega0::Float64=1.0)

    if e<TOLECC
        # select a smaller de for these, since we know we are expanding?


        # drop into circular frequency expansion calculations:
        Ω1 = Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e)
        β  = βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e) # = Ω2/Ω1
        Ω2 = β*Ω1

        Ω1plusa = Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a+da,e)
        Ω1pluse = Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e+de)

        ∂Ω1∂a = (Ω1plusa - Ω1)/da
        ∂Ω1∂e = (Ω1pluse - Ω1)/de

        ∂Ω2∂a = (Ω1plusa*βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a+da,e) - Ω2)/da
        ∂Ω2∂e = (Ω1pluse*βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e+de) - Ω2)/de

        return Ω1,Ω2,∂Ω1∂a,∂Ω1∂e,∂Ω2∂a,∂Ω2∂e

    else
        # using theta calculations to compute frequencies: leans heavily on Theta from Ufunc.jl
        # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

        # currently using Simpson's 1/3 rule for integration: requires that NINT be even.
        # @IMPROVE: we could use a better(?) integration scheme

        rperi,rapo = rpra_from_ae(a,e)

        function u8func(u::Float64)
            # push integration forward on eight different quantities:
            # 1. Theta(u)
            # 2. Theta(u)/r^2(u)
            # 3. Theta(u)*vr(u)
            # 4. dTheta(u)/da
            # 5. dTheta(u)/da
            # 6. dTheta(u)/da/r(u)^2
            # 7. Theta(u)/r^3(u)
            # 8. dTheta(u)/de/r(u)^2

            th = ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)
            dthda,dthde = ThetaAEdade(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE,da=da,de=de)

            return (th,
                    th/(ruAE(u,a,e)^2),
                    th*(Q(ψ,dψ,d2ψ,u,rperi,rapo)),
                    dthda,
                    dthde,
                    dthda/(ruAE(u,a,e)^2),
                    th/(ruAE(u,a,e)^3),
                    dthde/(ruAE(u,a,e)^2))

        end

        accum1,accum2,accum3,accum4,accum5,accum6,accum7,accum8 = UnitarySimpsonIntegration(u8func,K_O=NINT)

        #return the values
        Ω1inv = (1/pi)*accum1
        Ω1    = 1/Ω1inv
        actionj   = (1/pi)*accum3

        Eval, Lval, ∂E∂a, ∂E∂e, ∂L∂a, ∂L∂e = dELFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)

        # be careful with Ω2 if near radial: use analytic relation
        if e>(1-TOLECC)
            Ω2 = 0.5*Ω1
        else
            Ω2 = Lval*accum2*(1/pi)*Ω1
        end

        # now do the partial derivatives:
        ∂Ω1∂a = -((Ω1^2)/pi)*accum4
        ∂Ω1∂e = -((Ω1^2)/pi)*accum5

        # get E,L derivatives
        beta = Ω2/Ω1
        dbetada   = (∂L∂a*beta/Lval) - (2/a)*beta + (Lval/pi)*accum6
        #(1/pi) * ( (∂L∂a*beta) - (2*Lval/a)*beta + Lval*accum6)
        ∂Ω2∂a = (dbetada + (beta/Ω1)*∂Ω1∂a)*Ω1

        dbetade   = (∂L∂e*beta/Lval) - (2/e)*beta + (2*a*Lval/(e*pi))*accum7 + (Lval/pi)*accum8
        # (1/pi) * ( (∂L∂e*beta) - (2*Lval/e)*beta + (2*Lval*a/e)*accum7 + Lval*accum8)

        ∂Ω2∂e = (dbetade + (beta/Ω1)*∂Ω1∂e)*Ω1

        # return values: no option for action right now, but maybe?
        return Ω1,Ω2,∂Ω1∂a,∂Ω1∂e,∂Ω2∂a,∂Ω2∂e

    end # switches for orbits who are too radial or circular

end


"""DHenonThetaFreqRatiosAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,[da,de,action=false,NINT=32,EDGE=0.03,TOLECC=0.001])

use the defined function Theta(u) to compute frequency integrals
AND DERIVATIVES

returning alpha,beta and derivatives w.r.t. (a,e)

"""
function DHenonThetaFreqRatiosAE(ψ::Function,
                                 dψ::Function,
                                 d2ψ::Function,
                                 d3ψ::Function,
                                 d4ψ::Function,
                                 a::Float64,
                                 e::Float64;
                                 da::Float64=1.0e-6,
                                 de::Float64=1.0e-6,
                                 NINT::Int64=32,
                                 EDGE::Float64=0.01,
                                 TOLECC::Float64=0.001,
                                 Omega0::Float64=1.0)

    if e<TOLECC
        # drop into circular frequency expansion calculations:
        α  = Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e)/Omega0
        β  = βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)

        ∂α∂a = (Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a+da,e)/Omega0 - α)/da
        ∂α∂e = (Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e+de)/Omega0 - α)/de

        ∂β∂a = (βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a+da,e) - β)/da
        ∂β∂e = (βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e+de) - β)/de

        return α,β,∂α∂a,∂α∂e,∂β∂a,∂β∂e

    else
        # using theta calculations to compute frequencies: leans heavily on Theta from Ufunc.jl
        # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

        # currently using Simpson's 1/3 rule for integration: requires that NINT be even.
        # @IMPROVE: we could use a better(?) integration scheme

        rperi,rapo = rpra_from_ae(a,e)

        function u8func(u::Float64)
            # push integration forward on eight different quantities:
            # 1. Theta(u)
            # 2. Theta(u)/r^2(u)
            # 3. Theta(u)*vr(u)
            # 4. dTheta(u)/da
            # 5. dTheta(u)/da
            # 6. dTheta(u)/da/r(u)^2
            # 7. Theta(u)/r^3(u)
            # 8. dTheta(u)/de/r(u)^2

            th = ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)
            dthda,dthde = ThetaAEdade(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE,da=da,de=de)

            return (th,
                    th/(ruAE(u,a,e)^2),
                    th*(Q(ψ,dψ,d2ψ,u,rperi,rapo)),
                    dthda,
                    dthde,
                    dthda/(ruAE(u,a,e)^2),
                    th/(ruAE(u,a,e)^3),
                    dthde/(ruAE(u,a,e)^2))

        end

        accum1,accum2,accum3,accum4,accum5,accum6,accum7,accum8 = UnitarySimpsonIntegration(u8func,K_O=NINT)

        #return the values
        Ω1inv = (1/pi)*accum1
        Ω1    = 1/Ω1inv
        actionj   = (1/pi)*accum3

        Eval, Lval, ∂E∂a, ∂E∂e, ∂L∂a, ∂L∂e = dELFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)

        # be careful with Ω2 if near radial: use analytic relation
        if e>(1-TOLECC)
            Ω2 = 0.5*Ω1
        else
            Ω2 = Lval*accum2*(1/pi)*Ω1
        end

        # now do the partial derivatives:
        ∂α∂a = -((Ω1^2)/pi)*accum4/Omega0
        ∂α∂e = -((Ω1^2)/pi)*accum5/Omega0

        # get E,L derivatives
        β = Ω2/Ω1
        ∂β∂a   = (∂L∂a*β/Lval) - (2/a)*β + (Lval/pi)*accum6
        ∂β∂e   = (∂L∂e*β/Lval) - (2/e)*β + (2*a*Lval/(e*pi))*accum7 + (Lval/pi)*accum8

        # return values: no option for action right now, but maybe?
        return Ω1/Omega0,Ω2/Ω1,∂α∂a,∂α∂e,∂β∂a,∂β∂e

    end # switches for orbits who are too radial or circular

end






"""HenonThetaFrequenciesRpRa(ψ,dψ,d2ψ,rp,ra[,action=false,NINT=32,EDGE=0.03,TOLECC=0.001])

use the defined function Theta(u) to compute frequency integrals
"""
function HenonThetaFrequenciesRpRa(ψ::Function,
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

            th = ThetaRpRa(ψ,dψ,d2ψ,u,rperi,rapo,EDGE=EDGE)

            return (th,
                    th/(ruRpRa(u,rperi,rapo)^2),
                    th*(Q(ψ,dψ,d2ψ,u,rperi,rapo)))

        end

        accum1,accum2,accum3 = UnitarySimpsonIntegration(u3func,K_O=NINT)

        #return the values
        Ω1inv = (1/pi)*accum1
        Ω1    = 1/Ω1inv
        actionj   = (1/pi)*accum3

        # be careful with Ω2 if near radial: use analytic relation
        if e>(1-TOLECC)
            Ω2 = 0.5*Ω1
        else
            Ω2 = LFromRpRa(ψ,dψ,d2ψ,rperi,rapo)*accum2*(1/pi)*Ω1
        end

        if action
            return Ω1,Ω2,actionj
        else
            return Ω1,Ω2
        end

    end # switches for orbits who are too radial or circular

end




"""the henon anomaly increment
"""
function henon_f(u::Float64)
    u*(3/2 - (u^2)/2)
end

"""the derivative of the henon anomaly increment
"""
function henon_df(u::Float64)
    1.5*(1.0 - u^(2))
end

"""the second derivative of the henon anomaly increment
"""
function henon_d2f(u::Float64)
    return -3u
end

"""the third derivative of the henon anomaly increment
"""
function henon_d3f(u::Float64)
    return -3.
end

"""the fourth derivative of the henon anomaly increment
"""
function henon_d4f(u::Float64)
    return 0.
end


"""henon_anomaly_frequencies

use the henon anomaly mapping to compute orbit frequencies

@IMPROVE, handle negative sqrts gracefully
@IMPROVE, add switch for action calculation
@IMPROVE, offload transformations to OrbitDefinitions.jl
@IMPROVE, disallow any frequency overshoots with known boundaries

"""
function HenonAnomalyFrequencies(ψ::Function,
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


    #Ω1 = integration_distance/(accum1*dt);
    Ω1 = (pi/2)*uNINT/(accum1)
    Ω2 = Ω1/(pi/2) * L * (sm/am) * accum2 / uNINT;

    # @IMPROVE
    # note that we could also return the actions if we wanted:
    Jr = accum0/((pi/2)*uNINT)
    #action2 = jj;

    # @IMPROVE: we may want to force never allowing an overshoot (i.e. Ω1 is capped at 1 no matter what)
    if action
        return Ω1,Ω2,Jr
    else
        return Ω1,Ω2
    end

end


"""ComputeFrequenciesHenonRpRa(ψ,dψ/dr,d²ψ/dr²,rp,ra[,TOLECC,verbose])
"""
function ComputeFrequenciesHenonRpRa(ψ::Function,dψ::Function,d2ψ::Function,
        rperi::Float64,rapo::Float64;action::Bool=false,TOLECC::Float64=0.01,verbose::Int64=0,NINT::Int64=32)

    # check the tolerance
    a,ecc = ae_from_rpra(rperi,rapo)

    out = ComputeFrequenciesHenonAE(ψ,dψ,d2ψ,a,ecc,action=action,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

    if action
        return out[1],out[2],out[3]
    else
        return out[1],out[2]
    end

end

"""compute_frequences_henon_ae(ψ,dψ/dr,d²ψ/dr²,rp,ra[,TOLECC,verbose])
"""
function ComputeFrequenciesHenonAE(ψ::Function,dψ::Function,d2ψ::Function,
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

    # if radial, we should exploit Ω1 = 2*Ω2 somehow.

    # go to the frequency calculation
    freqs = HenonAnomalyFrequencies(ψ,rapo,rperi,E,J,action=true,NINT=NINT)

    if action
        Ω1,Ω2,Jr = freqs
        return Ω1,Ω2,Jr
    else
        Ω1,Ω2 = freqs
        return Ω1,Ω2
    end

end
