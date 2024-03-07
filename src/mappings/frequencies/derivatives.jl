
########################################################################
#
# (a,e) ↦ (α,β) mapping derivatives
# These methods for derivatives are not use/usable
# @IMPROVE: adapt to v2.0
# @IMPROVE: make them accessible while allowing for other integration methods
#
########################################################################

"""
    DαβHenonΘAE(ψ,dψ,d2ψ,a,e,params)

use the defined function Θ(u) to compute frequency ratios integrals
AND DERIVATIVES
"""
function DαβHenonΘAE(model::Potential,
                     a::Float64,e::Float64,
                     params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

    Ω₀ = params.Ω₀
    tola, tole = params.TOLA, EccentricityTolerance(a,params.rc,params.TOLECC)
    
    if (a < tola) || (e < tole) || (e > 1.0-tole)
        # For edge cases: Derivation outside the integral
        # Function to differentiate
        fun(atemp::Float64,etemp::Float64) = αβHenonΘAE(model,atemp,etemp,params)
        # Perform differentiation
        floc, ∂f∂a, ∂f∂e = NumericalDerivativeAE(fun,a,e,params)
        # Recast results
        α, β = floc
        ∂α∂a, ∂β∂a = ∂f∂a
        ∂α∂e, ∂β∂e = ∂f∂e

        return α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e
    end

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

        Θ = ΘAE(model,u,a,e,params)
        ∂Θ∂a, ∂Θ∂e = dΘAE(model,u,a,e,params)

        r = ru(u,a,e)

        return (Θ,
                Θ/(r^2),
                ∂Θ∂a,
                ∂Θ∂e,
                ∂Θ∂a/(r^2),
                (∂Θ∂e - 2.0*a*henonf(u)*Θ/r)/r^2)
    end

    accum1,accum2,accum3,accum4,accum5,accum6 = UnitarySimpsonIntegration(u6func,params.NINT)

    _, Lval, _, ∂L∂a, _, ∂L∂e = dELFromAE(model,a,e,params)

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
end

"""
    DFrequenciesHenonΘAE(ψ,dψ,d2ψ,a,e,params)

use the defined function Θ(u) to compute frequency integrals
AND DERIVATIVES
"""
function DFrequenciesHenonΘAE(model::Potential,
                              a::Float64,e::Float64,
                              params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

    Ω₀ = params.Ω₀
    α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = DαβHenonΘAE(model,a,e,params)

    Ω1, Ω2          = FrequenciesFromαβ(α,β,Ω₀)
    ∂Ω1∂a, ∂Ω2∂a    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂a,∂β∂a,Ω₀)
    ∂Ω1∂e, ∂Ω2∂e    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂e,∂β∂e,Ω₀)

    return Ω1, Ω2, ∂Ω1∂a, ∂Ω2∂a, ∂Ω1∂e, ∂Ω2∂e
end