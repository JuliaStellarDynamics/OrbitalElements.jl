

"""
HenonJFromAE(ψ,dψ,a,e,params)

compute the radial action alone
"""
function HenonJFromAE(ψ::F0,dψ::F1,
                      a::Float64,e::Float64,
                      params::OrbitalParameters=OrbitalParameters())::Float64 where {F0  <: Function, F1 <: Function}

    u1func(u::Float64)::Float64 = drdu(u,a,e)*Vrad(ψ,dψ,u,a,e,params)

    return (1/pi)*UnitarySimpsonIntegration(u1func,params.NINT)
end

"""
    αβHenonΘAE(ψ,dψ,d2ψ,a,e,params)

use the defined function Θ(u) to compute frequency ratios integrals.
"""
function αβHenonΘAE(ψ::F0,dψ::F1,d2ψ::F2,
                    a::Float64,e::Float64,
                    params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    Ω₀ = params.Ω₀
    tola, tole = params.TOLA, EccentricityTolerance(a,params.rc,params.TOLECC)
    if (a < tola)
        # 2nd order interpolation
        # between center (a=0) and value at a=tola and a=2*tola
        # Center:
        a0 = 0.
        α0, β0 = Ω1circular(dψ,d2ψ,a0), 0.5
        # e=tole:
        atola = tola
        αtola, βtola = αβHenonΘAE(ψ,dψ,d2ψ,atola,e,params)
        # e=2*tole:
        a2tola = 2*tola
        α2tola, β2tola = αβHenonΘAE(ψ,dψ,d2ψ,a2tola,e,params)

        # Interpolation
        α = Interpolation2ndOrder(a,a0,α0,atola,αtola,a2tola,α2tola)
        β = Interpolation2ndOrder(a,a0,β0,atola,βtola,a2tola,β2tola)

        return α, β
    elseif (e == 0.)
        Ω1, Ω2 = Ω1circular(dψ,d2ψ,a), Ω2circular(dψ,d2ψ,a)

        return αβFromFrequencies(Ω1,Ω2,Ω₀)
    elseif (0. < e < tole)
        # 2nd order interpolation
        # between circular (e=0) and value at e=tole and e=2*tole
        # Circular:
        ecirc = 0.
        αcirc, βcirc = αβHenonΘAE(ψ,dψ,d2ψ,a,ecirc,params)
        # e=tole:
        etole = tole
        αtole, βtole = αβHenonΘAE(ψ,dψ,d2ψ,a,etole,params)
        # e=2*tole:
        e2tole = 2*tole
        α2tole, β2tole = αβHenonΘAE(ψ,dψ,d2ψ,a,e2tole,params)

        # Interpolation
        α = Interpolation2ndOrder(e,ecirc,αcirc,etole,αtole,e2tole,α2tole)
        β = Interpolation2ndOrder(e,ecirc,βcirc,etole,βtole,e2tole,β2tole)

        return α, β
    elseif ((1.0-tole) < e < 1.)
        # 2nd order interpolation
        # between radial (e=1.) and value at e=1-tole and e=1-2*tole
        # Radial:
        erad = 1.
        αrad, βrad = αβHenonΘAE(ψ,dψ,d2ψ,a,erad,params)
        # e=1-tole:
        etole = 1.0-tole
        αtole, βtole = αβHenonΘAE(ψ,dψ,d2ψ,a,etole,params)
        # e=1-2*tole:
        e2tole = 1.0-2*tole
        α2tole, β2tole = αβHenonΘAE(ψ,dψ,d2ψ,a,e2tole,params)

        # Interpolation
        α = Interpolation2ndOrder(e,e2tole,α2tole,etole,αtole,erad,αrad)
        β = Interpolation2ndOrder(e,e2tole,β2tole,etole,βtole,erad,βrad)

        return α, β
    else
        # using Θ calculations to compute frequencies: leans heavily on Θ from Ufunc.jl
        # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

        # currently using Simpson's 1/3 rule for integration: requires that NINT be even.
        # Is it true ?? Why ?
        # @IMPROVE: we could use a better integration scheme

        # this doesn't throw any allocations, so don't worry about that!
        function u2func(u::Float64)::Tuple{Float64,Float64}
            # push integration forward on two different quantities: Θ(u),Θ(u)/r^2(u)
            Θ = ΘAE(ψ,dψ,d2ψ,u,a,e,params)

            return Θ, Θ/(ru(u,a,e)^2)
        end

        accum1, accum2 = UnitarySimpsonIntegration(u2func,params.NINT)

        #return the values
        invα = (Ω₀/pi)*accum1
        α    = 1.0/invα
        β = (e == 1.) ? 0.5 : (LFromAE(ψ,dψ,a,e,params)/pi)*accum2

        return α, β
    end # switches for orbits who are too radial or circular
end

"""
    DαβHenonΘAE(ψ,dψ,d2ψ,a,e,params)

use the defined function Θ(u) to compute frequency ratios integrals
AND DERIVATIVES
"""
function DαβHenonΘAE(ψ::F0,dψ::F1,d2ψ::F2,
                     a::Float64,e::Float64,
                     params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    Ω₀ = params.Ω₀
    tola, tole = params.TOLA, EccentricityTolerance(a,params.rc,params.TOLECC)
    # Numerical derivative points
    ap, da, ep, de = NumDerivPoints(a,e,params.da,params.de,tola,tole)
    
    if (a < tola) || (e < tole) || (e > 1.0-tole)
        # For edge cases: Derivation outside the integral
        # Current point
        α, β = αβHenonΘAE(ψ,dψ,d2ψ,a,e,params)

        # For a derivatives
        αap, βap = αβHenonΘAE(ψ,dψ,d2ψ,ap,e,params)
        
        # For e derivatives
        αep, βep = αβHenonΘAE(ψ,dψ,d2ψ,a,ep,params)

        ∂α∂a = (αap-α)/da
        ∂β∂a = (βap-β)/da

        ∂α∂e = (αep-α)/de
        ∂β∂e = (βep-β)/de

        return α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e
    else
        # Derivation inside the integral
        # using Θ calculations to compute frequencies: leans heavily on Θ from Ufunc.jl
        # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

        # WARNING !! Strong assumption:
        # r(u) = a(1+ef(u))
        function u6func(u::Float64)
            # push integration forward on eight different quantities:
            # 1. Θ                          → α
            # 2. Θ/r^2                      → β
            # 3. ∂Θ/∂a                      → ∂α∂a
            # 4. ∂Θ/∂e                      → ∂α∂e
            # 5. ∂Θ/∂a/r^2                  → ∂β∂a
            # 6. (∂Θ/∂e - 2af(u)Θ/r )/r^2   → ∂β∂e

            Θ = ΘAE(ψ,dψ,d2ψ,u,a,e,params)
            ∂Θ∂a, ∂Θ∂e = dΘAE(ψ,dψ,d2ψ,u,a,e,params)

            r = ru(u,a,e)

            return (Θ,
                    Θ/(r^2),
                    ∂Θ∂a,
                    ∂Θ∂e,
                    ∂Θ∂a/(r^2),
                    (∂Θ∂e - 2.0*a*henonf(u)*Θ/r)/r^2)
        end

        accum1,accum2,accum3,accum4,accum5,accum6 = UnitarySimpsonIntegration(u6func,params.NINT)

        _, Lval, _, ∂L∂a, _, ∂L∂e = dELFromAE(ψ,dψ,a,e,params)

        # α
        invα = (Ω₀/pi)*accum1
        α    = 1.0/invα

        β = (Lval/pi)*accum2

        # ∂α
        ∂α∂a = -(α^2)*(Ω₀/pi)*accum3
        ∂α∂e = -(α^2)*(Ω₀/pi)*accum4

        # ∂β
        βoverL = accum2/pi
        ∂β∂a = ∂L∂a*βoverL - 2.0*β/a + (Lval/pi)*accum5
        ∂β∂e = ∂L∂e*βoverL + (Lval/pi)*accum6

        return α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e
    end # switches for orbits who are too radial or circular
end

"""
    DFrequenciesHenonΘAE(ψ,dψ,d2ψ,a,e,params)

use the defined function Θ(u) to compute frequency integrals
AND DERIVATIVES
"""
function DFrequenciesHenonΘAE(ψ::F0,dψ::F1,d2ψ::F2,
                              a::Float64,e::Float64,
                              params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    Ω₀ = params.Ω₀
    α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = DαβHenonΘAE(ψ,dψ,d2ψ,a,e,params)

    Ω1, Ω2          = FrequenciesFromαβ(α,β,Ω₀)
    ∂Ω1∂a, ∂Ω2∂a    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂a,∂β∂a,Ω₀)
    ∂Ω1∂e, ∂Ω2∂e    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂e,∂β∂e,Ω₀)

    return Ω1, Ω2, ∂Ω1∂a, ∂Ω2∂a, ∂Ω1∂e, ∂Ω2∂e
end