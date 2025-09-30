
########################################################################
#
# (a,e) ↦ (α,β) mapping derivatives
# These methods for derivatives are not use/usable
# @IMPROVE: adapt to v2.0
# @IMPROVE: make them accessible while allowing for other integration methods
# @IMPROVE: Blue style
########################################################################

"""
    αβ_from_ae_internal_derivatives(a, e, model, params)

use the defined function _Θ(w) to compute frequency ratios integrals
AND DERIVATIVES
"""
function αβ_from_ae_internal_derivatives(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}
    tola, tole = params.TOLA, _eccentricity_tolerance(a/params.rc, params.TOLECC)
    if (a < tola) || (e < tole) || (e > 1.0-tole)
        # For edge cases: Derivation outside the integral
        # Function to differentiate
        fun(atmp::Float64,etmp::Float64) = αβ_from_AE(atmp, etmp, model, params)
        # Perform differentiation
        floc, ∂f∂a, ∂f∂e = _derivatives_ae(fun, a, e, params)
        # Recast results
        α, β = floc
        ∂α∂a, ∂β∂a = ∂f∂a
        ∂α∂e, ∂β∂e = ∂f∂e

        return α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e
    end
    # Derivation inside the integral
    # using Θ calculations to compute frequencies: leans heavily on _Θ from 
    # radial_velocity.jl
    # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi

    # @WARNING !! Strong assumption:
    # r(w) = a(1+e * _henonf(w))
    function w6func(w::Float64)
        # push integration forward on eight different quantities:
        # 1. Θ                          → α
        # 2. Θ/r^2                      → β
        # 3. ∂Θ/∂a                      → ∂α∂a
        # 4. ∂Θ/∂e                      → ∂α∂e
        # 5. ∂Θ/∂a/r^2                  → ∂β∂a
        # 6. (∂Θ/∂e - 2af(w)Θ/r )/r^2   → ∂β∂e

        integrand = _Θ(w, a, e, model, params)
        ∂Θ∂a, ∂Θ∂e = _Θ_derivatives_ae(w, a, e, model, params)

        r = radius_from_anomaly(w, a, e, model, params)

        return (
            integrand,
            integrand / (r^2),
            ∂Θ∂a,
            ∂Θ∂e,
            ∂Θ∂a / (r^2),
            (∂Θ∂e - 2a * _henonf(w) * Θ / r) / (r^2)
        )
    end

    accum1, accum2, accum3, accum4, accum5, accum6 = _integrate_simpson(w6func, params.NINT)
    _, Lval, _, ∂L∂a, _, ∂L∂e = EL_from_ae_derivatives(a, e, model, params)
    # α
    Ω0 = frequency_scale(model)
    invα = (Ω0 / pi) * accum1
    α    = 1.0/invα
    β = (Lval/pi)*accum2
    # ∂α
    ∂α∂a = -(α^2) * (Ω0 / pi) * accum3
    ∂α∂e = -(α^2) * (Ω0 / pi) * accum4
    # ∂β
    βoverL = accum2 / pi
    ∂β∂a = ∂L∂a * βoverL - 2β / a + (Lval / pi) * accum5
    ∂β∂e = ∂L∂e * βoverL + (Lval / pi) * accum6

    return α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e
end