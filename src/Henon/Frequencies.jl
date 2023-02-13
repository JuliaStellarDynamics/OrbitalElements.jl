

"""
    HenonJFromAE(ψ,dψ,d2ψ,d3ψ,a,e,NINT,TOLECC)

compute the radial action alone
"""
function HenonJFromAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                      a::Float64,e::Float64,
                      NINT::Int64,TOLECC::Float64)::Float64 where {F0  <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    u1func(u::Float64)::Float64 = drdu(u,a,e)*Vrad(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLECC)

    return (1/pi)*UnitarySimpsonIntegration(u1func,NINT)
end


"""
    αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)

use the defined function Θ(u) to compute frequency integrals
"""
function αβHenonΘAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                    a::Float64,e::Float64,
                    TOLA::Float64,TOLECC::Float64,
                    NINT::Int64,EDGE::Float64,Ω₀::Float64=1.0)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    tole = EccentricityTolerance(a,TOLA,TOLECC)
    if (e <= tole)
        # drop into circular frequency expansion calculations:
        α = αcircular(dψ,d2ψ,d3ψ,d4ψ,a,e,Ω₀)
        β  = βcircular(dψ,d2ψ,d3ψ,d4ψ,a,e) # = Ω2/Ω1
        #Ω2 = β*Ω1

        return α,β

    elseif (1.0-tole) < e < 1.
        # No taylor expansion for this radial orbits case
        # β fixed at 0.5 quite poor approximation (especially for derivatives)
        # → linear interpolation between 1/2 (e=1) and value at e=1.0-TOLECC
        # βtole (quasi-radial)
        etole = 1.0 - tole
        αtole, βtole = αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,etole,TOLA,TOLECC,NINT,EDGE,Ω₀)
        # For α value in e = 1. depends on the potential (non trivial dependence).
        # The value computed through Θ integration in e = 1. is not absurd
        # Taking linear interpolation between 1-tole and 1 prevents
        # huge error on the α(e) curve's slope (hence better for the derivatives)
        erad = 1.0
        αrad, βrad = αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,erad,TOLA,TOLECC,NINT,EDGE,Ω₀)

        e2tole = 1.0 - 2*tole
        α2tole, β2tole = αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e2tole,TOLA,TOLECC,NINT,EDGE,Ω₀)

        # Interpolation
        α = Interpolation2ndOrder(e,e2tole,α2tole,etole,αtole,erad,αrad)
        β = Interpolation2ndOrder(e,e2tole,β2tole,etole,βtole,erad,βrad)

        return α, β

    else

        # using Θ calculations to compute frequencies: leans heavily on Θ from Ufunc.jl
        # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

        # currently using Simpson's 1/3 rule for integration: requires that NINT be even.
        # @IMPROVE: we could use a better integration scheme

        # this doesn't throw any allocations, so don't worry about that!
        function u2func(u::Float64)::Tuple{Float64,Float64}
            # push integration forward on three different quantities: Θ(u),Θ(u)/r^2(u)
            Θ = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC,EDGE)

            return Θ, Θ/(ru(u,a,e)^2)
        end

        accum1, accum2 = UnitarySimpsonIntegration(u2func,NINT)

        #return the values
        #return the values
        invα = (Ω₀/pi)*accum1
        α    = 1.0/invα
        β = (e == 1.) ? 0.5 : (LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC)/pi)*accum2

        return α, β

    end # switches for orbits who are too radial or circular

end


"""
    DαβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)

use the defined function Θ(u) to compute frequency ratios integrals
AND DERIVATIVES
"""
function DαβHenonΘAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                     a::Float64,e::Float64,
                     da::Float64=1.0e-6,de::Float64=1.0e-6,
                     TOLA::Float64=0.001,TOLECC::Float64=0.001,
                     NINT::Int64=32,EDGE::Float64=0.01,
                     Ω₀::Float64=1.0)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    tole = EccentricityTolerance(a,TOLA,TOLECC)
    # Numerical derivative points
    ap, da, ep, de = NumDerivPoints(a,e,da,de,TOLA,tole)

    # if nearly circular, use the epicyclic approximation
    if e <= tole
        # drop into circular frequency expansion calculations:
        # Derivation outside the integral !
        α = αcircular(dψ,d2ψ,d3ψ,d4ψ,a,e,Ω₀)
        β  = βcircular(dψ,d2ψ,d3ψ,d4ψ,a,e)

        # For a derivatives
        αap = αcircular(dψ,d2ψ,d3ψ,d4ψ,ap,e,Ω₀)
        βap  = βcircular(dψ,d2ψ,d3ψ,d4ψ,ap,e)

        # For e derivatives
        αep = αcircular(dψ,d2ψ,d3ψ,d4ψ,a,ep,Ω₀)
        βep  = βcircular(dψ,d2ψ,d3ψ,d4ψ,a,ep)

        ∂α∂a = (αap-α)/da
        ∂β∂a = (βap-β)/da

        ∂α∂e = (αep-α)/de
        ∂β∂e = (βep-β)/de

        return α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e
    elseif e>(1.0-tole)
        # Derivation outside the integral
        α  , β   = αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLECC,NINT,EDGE,Ω₀)

        # For a derivatives
        αap, βap = αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,ap,e,TOLECC,NINT,EDGE,Ω₀)

        # For e derivatives
        αep, βep = αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,ep,TOLECC,NINT,EDGE,Ω₀)

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

            Θ = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC,EDGE)
            ∂Θ∂a, ∂Θ∂e = dΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,da,de,TOLA,TOLECC,EDGE)

            r = ru(u,a,e)

            return (Θ,
                    Θ/(r^2),
                    ∂Θ∂a,
                    ∂Θ∂e,
                    ∂Θ∂a/(r^2),
                    (∂Θ∂e - 2.0*a*henon_f(u)*Θ/r)/r^2)
        end

        accum1,accum2,accum3,accum4,accum5,accum6 = UnitarySimpsonIntegration(u6func,NINT)

        _, Lval, _, ∂L∂a, _, ∂L∂e = dELFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC)

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
    DFrequenciesHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)

use the defined function Θ(u) to compute frequency integrals
AND DERIVATIVES
"""
function DFrequenciesHenonΘAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              a::Float64,e::Float64,
                              da::Float64=1.0e-6,de::Float64=1.0e-6,
                              TOLA::Float64=0.001,TOLECC::Float64=0.001,
                              NINT::Int64=32,EDGE::Float64=0.01,
                              Ω₀::Float64=1.0)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = DαβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)

    Ω1, Ω2          = FrequenciesFromαβ(α,β,Ω₀)
    ∂Ω1∂a, ∂Ω2∂a    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂a,∂β∂a,Ω₀)
    ∂Ω1∂e, ∂Ω2∂e    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂e,∂β∂e,Ω₀)

    return Ω1, Ω2, ∂Ω1∂a, ∂Ω2∂a, ∂Ω1∂e, ∂Ω2∂e
end
