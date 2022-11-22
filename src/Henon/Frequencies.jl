

"""HenonΘFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e,[,action=false,NINT=32,EDGE=0.03,TOLECC=0.001])

use the defined function Θ(u) to compute frequency integrals

this is the fast version, without action computation, and more parameters specified
"""
@inline function HenonΘFrequenciesAE(ψ::Function,
                                 dψ::Function,
                                 d2ψ::Function,
                                 d3ψ::Function,
                                 d4ψ::Function,
                                 a::Float64,
                                 e::Float64,
                                 TOLECC::Float64,
                                 NINT::Int64,
                                 EDGE::Float64)::Tuple{Float64,Float64}

    if e<TOLECC

        #(VERBOSE > 1) && println("OrbitalElements.Henon.Frequencies.HenonΘFrequenciesAE: using circular approximation for a=$a, e=$e.")

        # drop into circular frequency expansion calculations:
        Ω1 = Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e)
        β  = βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e) # = Ω2/Ω1
        Ω2 = β*Ω1

        return Ω1,Ω2

    else

        # using Θ calculations to compute frequencies: leans heavily on Θ from Ufunc.jl
        # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

        # currently using Simpson's 1/3 rule for integration: requires that NINT be even.
        # @IMPROVE: we could use a better integration scheme

        # this doesn't throw any allocations, so don't worry about that!
        function u2func(u::Float64)::Tuple{Float64,Float64}
            # push integration forward on three different quantities: Θ(u),Θ(u)/r^2(u)

            #th = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLECC,EDGE=EDGE)
            th = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLECC,EDGE)

            return (th,
                    th/(ru(u,a,e)^2))

        end

        accum = UnitarySimpsonIntegration(u2func,NINT)

        #return the values
        Ω1inv = (1/pi)*accum[1]
        Ω1    = 1/Ω1inv

        # be careful with Ω2 if near radial: use analytic relation
        if e>(1.0-TOLECC)
            Ω2 = 0.5*Ω1
        else
            Ω2 = LFromAE(ψ,dψ,d2ψ,d3ψ,a,e)*accum[2]*(1/pi)*Ω1
        end

        return Ω1,Ω2

    end # switches for orbits who are too radial or circular

end

"""
for use when computing the radial action
"""
function HenonΘJFrequenciesAE(ψ::Function,
                                 dψ::Function,
                                 d2ψ::Function,
                                 d3ψ::Function,
                                 d4ψ::Function,
                                 a::Float64,
                                 e::Float64;
                                 action::Bool=false,
                                 NINT::Int64=32,
                                 EDGE::Float64=0.03,
                                 TOLECC::Float64=ELTOLECC)::Tuple{Float64,Float64,Float64}

    if e<TOLECC

        #(VERBOSE > 1) && println("OrbitalElements.Henon.Frequencies.HenonΘFrequenciesAE: using circular approximation for a=$a, e=$e.")

        # drop into circular frequency expansion calculations:
        Ω1 = Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e)
        β  = βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e) # = Ω2/Ω1
        Ω2 = β*Ω1

        if action
            u1func(u::Float64)::Float64 = drdu(u,a,e)*Vrad(ψ,dψ,d2ψ,d3ψ,u,a,e)
            accum = UnitarySimpsonIntegration(u1func,NINT)
            actionj   = (1/pi)*accum
            return Ω1,Ω2,actionj
        else
            return Ω1,Ω2
        end

    else

        # using Θ calculations to compute frequencies: leans heavily on Θ from Ufunc.jl
        # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

        # currently using Simpson's 1/3 rule for integration: requires that NINT be even.
        # @IMPROVE: we could use a better integration scheme

        function u3func(u::Float64)
            # push integration forward on three different quantities: Θ(u),Θ(u)/r^2(u),Θ(u)*vr(u)

            #th = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLECC,EDGE=EDGE)
            th = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLECC,EDGE)

            return (th,
                    th/(ru(u,a,e)^2),
                    drdu(u,a,e)*Vrad(ψ,dψ,d2ψ,d3ψ,u,a,e))

        end

        accum = UnitarySimpsonIntegration(u3func,NINT)

        #return the values
        Ω1inv = (1/pi)*accum[1]
        Ω1    = 1/Ω1inv
        actionj   = (1/pi)*accum[3]

        # be careful with Ω2 if near radial: use analytic relation
        if e > (1.0-TOLECC)
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

"""DHenonΘFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,[da,de,action=false,NINT=32,EDGE=0.03,TOLECC=0.001])

use the defined function Θ(u) to compute frequency integrals
AND DERIVATIVES
"""
function DHenonΘFrequenciesAE(ψ::Function,
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

    # if nearly circular, use the epicyclic approximation
    if e<TOLECC

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
        # using Θ calculations to compute frequencies: leans heavily on Θ from Ufunc.jl
        # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

        # currently using Simpson's 1/3 rule for integration: requires that NINT be even.
        # @IMPROVE: we could use a better(?) integration scheme

        function u8func(u::Float64)
            # push integration forward on eight different quantities:
            # 1. Θ(u)
            # 2. Θ(u)/r^2(u)
            # 3. Θ(u)*(vr(u))^2
            # 4. dΘ(u)/da
            # 5. dΘ(u)/da
            # 6. dΘ(u)/da/r(u)^2
            # 7. Θ(u)/r^3(u)
            # 8. dΘ(u)/de/r(u)^2

            #th = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE)
            th = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,ELTOLECC,EDGE)
            #dthda,dthde = ΘAEdade(ψ,dψ,d2ψ,d3ψ,u,a,e,EDGE=EDGE,da=da,de=de)
            dthda,dthde = ΘAEdade(ψ,dψ,d2ψ,d3ψ,u,a,e,ELTOLECC,EDGE,da,de)

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

        accum1,accum2,accum3,accum4,accum5,accum6,accum7,accum8 = UnitarySimpsonIntegration(u8func,NINT)

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
        β = Ω2/Ω1
        dβda   = (∂L∂a*β/Lval) - (2/a)*β + (Lval/pi)*accum6
        ∂Ω2∂a = (dβda + (β/Ω1)*∂Ω1∂a)*Ω1

        dβde   = (∂L∂e*β/Lval) - (2/e)*β + (2*a*Lval/(e*pi))*accum7 + (Lval/pi)*accum8
        ∂Ω2∂e = (dβde + (β/Ω1)*∂Ω1∂e)*Ω1

        # return values: no option for action right now, but maybe?
        return Ω1,Ω2,∂Ω1∂a,∂Ω1∂e,∂Ω2∂a,∂Ω2∂e

    end # switches for orbits who are too radial or circular
end


"""DHenonΘFreqRatiosAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,[da,de,action=false,NINT=32,EDGE=0.03,TOLECC=0.001,Ω₀=1.0])
returning α,β and derivatives w.r.t. (a,e) using DHenonΘFrequenciesAE
"""
function DHenonΘFreqRatiosAE(ψ::Function,
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
                             Ω₀::Float64=1.0)

    # Frenquency computatio (with derivatives)
    #Ω1,Ω2,∂Ω1∂a,∂Ω1∂e,∂Ω2∂a,∂Ω2∂e = DHenonΘFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;da=da,de=de,NINT=NINT,EDGE=EDGE,TOLECC=TOLECC)
    Ω1,Ω2,∂Ω1∂a,∂Ω1∂e,∂Ω2∂a,∂Ω2∂e = DHenonΘFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;da=da,de=de,NINT=NINT,EDGE=EDGE,TOLECC=TOLECC)

    # From (Ω1,Ω2) to (α,β)
    α = Ω1 / Ω₀
    ∂α∂a = ∂Ω1∂a / Ω₀
    ∂α∂e = ∂Ω1∂e / Ω₀

    β = Ω2/Ω1
    ∂β∂a   = 1.0 / Ω1 * (∂Ω2∂a - β*∂Ω1∂a)
    ∂β∂e   = 1.0 / Ω1 * (∂Ω2∂e - β*∂Ω1∂e)

    # return values: no option for action right now, but maybe?
    return α,β,∂α∂a,∂α∂e,∂β∂a,∂β∂e
end

function HenonJFromAE(ψ::Function,
                      dψ::Function,
                      d2ψ::Function,
                      d3ψ::Function,
                      a::Float64,
                      e::Float64;
                      NINT::Int64=32)

    u1func(u::Float64)::Float64 = drdu(u,a,e)*Vrad(ψ,dψ,d2ψ,d3ψ,u,a,e)

    return (1/pi)*UnitarySimpsonIntegration(u1func,NINT)
end
