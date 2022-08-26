

"""HenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e,[,action=false,NINT=32,EDGE=0.03,TOLECC=0.001])

use the defined function Theta(u) to compute frequency integrals
"""
function HenonThetaFrequenciesAE(ψ::Function,
                                 dψ::Function,
                                 d2ψ::Function,
                                 d3ψ::Function,
                                 d4ψ::Function,
                                 a::Float64,
                                 e::Float64;
                                 action::Bool=false,
                                 NINT::Int64=32,
                                 EDGE::Float64=0.03,
                                 TOLECC::Float64=0.001,
                                 verbose::Int64=0)

    if e<TOLECC

        if verbose > 1
            println("Henon/Frequencies.jl:HenonThetaFrequenciesAE: using circular approximation.")
        end

        # drop into circular frequency expansion calculations:
        Ω1 = Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e)
        β  = βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e) # = Ω2/Ω1
        Ω2 = β*Ω1

        if action
            u1func(u::Float64)::Float64 = drdu(u,a,e)*Vrad(ψ,dψ,d2ψ,d3ψ,u,a,e)
            accum = UnitarySimpsonIntegration(u1func,K_O=NINT)
            actionj   = (1/pi)*accum
            return Ω1,Ω2,actionj
        else
            return Ω1,Ω2
        end

    else

        # using theta calculations to compute frequencies: leans heavily on Theta from Ufunc.jl
        # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

        # currently using Simpson's 1/3 rule for integration: requires that NINT be even.
        # @IMPROVE: we could use a better integration scheme

        function u3func(u::Float64)
            # push integration forward on three different quantities: Theta(u),Theta(u)/r^2(u),Theta(u)*vr(u)

            th = ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)

            return (th,
                    th/(ru(u,a,e)^2),
                    drdu(u,a,e)*Vrad(ψ,dψ,d2ψ,d3ψ,u,a,e))

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
                                  TOLECC::Float64=0.001)

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

        function u8func(u::Float64)
            # push integration forward on eight different quantities:
            # 1. Theta(u)
            # 2. Theta(u)/r^2(u)
            # 3. Theta(u)*(vr(u))^2
            # 4. dTheta(u)/da
            # 5. dTheta(u)/da
            # 6. dTheta(u)/da/r(u)^2
            # 7. Theta(u)/r^3(u)
            # 8. dTheta(u)/de/r(u)^2

            th = ThetaAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)
            dthda,dthde = ThetaAEdade(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE,da=da,de=de)

            r = ru(u,a,e)

            return (th,
                    th/(r^2),
                    drdu(u,a,e)*Vrad(ψ,dψ,d2ψ,d3ψ,u,a,e),
                    dthda,
                    dthde,
                    dthda/(r^2),
                    th/(r^3),
                    dthde/(r^2))

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


"""DHenonThetaFreqRatiosAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,[da,de,action=false,NINT=32,EDGE=0.03,TOLECC=0.001,Ω0=1.0])
returning alpha,beta and derivatives w.r.t. (a,e) using DHenonThetaFrequenciesAE
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
                                 Ω0::Float64=1.0)

    # Frenquency computatio (with derivatives)
    Ω1,Ω2,∂Ω1∂a,∂Ω1∂e,∂Ω2∂a,∂Ω2∂e = DHenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;da=da,de=de,NINT=NINT,EDGE=EDGE,TOLECC=TOLECC)

    # From (Ω1,Ω2) to (α,β)
    α = Ω1 / Ω0
    ∂α∂a = ∂Ω1∂a / Ω0
    ∂α∂e = ∂Ω1∂e / Ω0

    β = Ω2/Ω1
    ∂β∂a   = 1.0 / Ω1 * (∂Ω2∂a - β*∂Ω1∂a)
    ∂β∂e   = 1.0 / Ω1 * (∂Ω2∂e - β*∂Ω1∂e)

    # return values: no option for action right now, but maybe?
    return α,β,∂α∂a,∂α∂e,∂β∂a,∂β∂e
end