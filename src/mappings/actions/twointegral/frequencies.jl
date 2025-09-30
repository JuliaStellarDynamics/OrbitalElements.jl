"""
Frequency mapping definitions

@IMPROVE: Hénon anomaly and integrator (RK4) are hard-
coded. Should smartly incorparate them as parameters in the OrbitalParameters structure.
@IMPROVE: Would not work for cuspy potential due to the values at ``a=0``, but should use 
the same function in regular places ! Need to handle edges differently in this case. Be very
careful about infinite loop between edge and regular functions when implementing a solution.
"""

########################################################################
#
# (α,β) ↔ (Ω1,Ω2) mappings
#
########################################################################
"""
    αβ_from_frequencies(Ω1, Ω2, model)

frequencies ratios from frequencies
"""
function αβ_from_frequencies(Ω1::Float64, Ω2::Float64, Ω0::Float64)::Tuple{Float64,Float64}
    return Ω1 / Ω0, Ω2 / Ω1
end
function αβ_from_frequencies(
    Ω1::Float64,
    Ω2::Float64,
    model::TwoIntegralPotential
)::Tuple{Float64,Float64}
    return αβ_from_frequencies(Ω1, Ω2, frequency_scale(model))
end


"""
    frequencies_from_αβ(α, β, model)

frequencies from frequency ratios
"""
function frequencies_from_αβ(α::Float64, β::Float64, Ω0::Float64)::Tuple{Float64,Float64}
    return Ω0 * α, Ω0 * α * β
end
function frequencies_from_αβ(
    Ω1::Float64,
    Ω2::Float64,
    model::TwoIntegralPotential
)::Tuple{Float64,Float64}
    return frequencies_from_αβ(Ω1, Ω2, frequency_scale(model))
end

########################################################################
#
# (α,β) ↔ (Ω1,Ω2) mappings: derivatives and jacobian
#
########################################################################
"""
    αβ_from_frequencies_derivatives(Ω1, Ω2, model)

same as `αβ_from_frequencies` plus mapping derivatives.
"""
function αβ_from_frequencies_derivatives(
    Ω1::Float64,
    Ω2::Float64,
    model::TwoIntegralPotential
)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

    α, β = αβ_from_frequencies(Ω1, Ω2, model)
    ∂α∂Ω1, ∂β∂Ω1 = 1 / frequency_scale(model), -Ω2 / Ω1^2
    ∂α∂Ω2, ∂β∂Ω2 = 0.0, 1 / Ω1
    return α, β, ∂α∂Ω1, ∂β∂Ω1, ∂α∂Ω2, ∂β∂Ω2
end

"""
    frequencies_from_αβ_derivatives(α, β, model)

same as `frequencies_from_αβ` plus mapping derivatives.
"""
function frequencies_from_αβ_derivatives(
    α::Float64,
    β::Float64,
    model::TwoIntegralPotential
)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

    Ω0 = frequency_scale(model)
    Ω1, Ω2 = frequencies_from_αβ(α, β, Ω0)
    ∂Ω1∂α, ∂Ω2∂α = Ω0 .* (1, β)
    ∂Ω1∂β, ∂Ω2∂β = Ω0 .* (0, α)
    return Ω1, Ω2, ∂Ω1∂α, ∂Ω2∂α, ∂Ω1∂β, ∂Ω2∂β
end

"""
    frequencies_to_αβ_jacobian(Ω1, Ω2, model)

jacobian of the (Ω1,Ω2) ↦ (α,β) mapping, i.e. |∂(α,β)/∂(Ω1,Ω2)|.
"""
function frequencies_to_αβ_jacobian(
    Ω1::Float64,
    Ω2::Float64,
    model::TwoIntegralPotential
)
    return abs(1 / (frequency_scale(model) * Ω1))
end

"""
    αβ_to_frequencies_jacobian(α, β, model)

jacobian of the (α,β) ↦ (Ω1,Ω2) mapping, i.e. |∂(Ω1,Ω2)/∂(α,β)|.
"""
function αβ_to_frequencies_jacobian(
    α::Float64,
    β::Float64,
    model::TwoIntegralPotential
)
    return abs(frequency_scale(model)^2 * α)
end

########################################################################
#
# (a,e) ↦ frequencies mapping 
#
#######################################################################
include("frequencies/anomaly.jl")
include("frequencies/radial_velocity.jl")
include("frequencies/circular.jl")

"""
    αβ_from_ae(a, e, model[, params])

mapping from semimajor axis and eccentricity to frequency ratios ``(\\alpha,\\beta)``.

@IMPROVE: Hénon anomaly and integrator (RK4) are hard-
coded. Should smartly incorparate them as parameters in the OrbitalParameters structure.
@IMPROVE: Would not work for cuspy potential due to the values at ``a=0``, but should use 
the same function in regular places ! Need to handle edges differently in this case. Be very
careful about infinite loop between edge and regular functions when implementing a solution.
"""
function αβ_from_ae(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # Edge cases
    if a == 0 || e == 0 # Inner or circular frequencies
        return _αcircular(a, model, params), _βcircular(a, model, params)
    end
    # Handling (a,e)-domain edges through interpolation
    # IMPORTANT : has to be before the generic computation
    fun(atmp::Float64, etmp::Float64) = αβ_from_ae(atmp, etmp, model, params)
    res = _interpolate_edges_ae(fun, a, e, params)
    if !isnothing(res)
        return res
    end
    # Generic computations
    # 1/α = Ω₀ / π * \int_{-1}^1 dw Θ(w); β = L / π * \int_{-1}^1 dw Θ(w) / [r(w)]^2
    # using Θ calculations to compute frequencies: leans heavily on _Θ from 
    # frequencies/anomaly.jl
    # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi
    # currently using Simpson's 1/3 rule only (hard-coded)
    # @IMPROVE: free the integration scheme, i.e., make it a parameter
    function invαβ_integrands(w::Float64)::Tuple{Float64,Float64}
        # push integration forward on two different quantities: Θ(w), Θ(w)/r^2(w)
        integrand = _Θ(w, a, e, model, params)
        r = radius_from_anomaly(w, a, e, model, params)
        return integrand, integrand / r^2
    end
    accum1, accum2 = _integrate_simpson(invαβ_integrands, params.NINT)
    _, L = EL_from_ae(a, e, model, params)
    α = pi / (frequency_scale(model) * accum1)
    β = e == 1 ? 1/2 : (L / pi) * accum2
    return α, β
end

"""
    frequencies_from_ae(a,e,model,params)

mapping from semimajor axis and eccentricity to frequencies ``(\\Omega_1,\\Omega_2)``.
"""
function frequencies_from_ae(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    α, β = αβ_from_ae(a, e, model, params)
    return frequencies_from_αβ(α, β, model)
end


########################################################################
#
# (a,e) ↦ frequencies mapping: derivatives
#
########################################################################
"""
    αβ_from_ae_derivatives(a, e, model[, params])

same as `αβ_from_ae` plus mapping derivatives.

@IMPROVE: right now the edge handling procedure is hard-coded. Should incorparate it in the
OrbitalParameters structure (not only the tolerances)
@IMPROVE: right now the derivative procedure is hard-coded. Should incorparate it 
in the OrbitalParameters structure (not only the steps)
"""
function αβ_from_ae_derivatives(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}
    # Function to differentiate
    fun(atmp::Float64,etmp::Float64) = αβ_from_ae(atmp, etmp, model, params)
    # Perform differentiation
    floc, ∂f∂a, ∂f∂e = _derivatives_ae(fun, a, e, params)
    # Recast results
    α, β = floc
    ∂α∂a, ∂β∂a = ∂f∂a
    ∂α∂e, ∂β∂e = ∂f∂e

    return α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e
end

"""
    frequencies_from_ae_derivatives(a, e, model[, params])
"""
function frequencies_from_ae_derivatives(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}
    α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = αβ_from_ae_derivatives(a, e, model, params)
    Ω1, Ω2, ∂Ω1∂α, ∂Ω2∂α, ∂Ω1∂β, ∂Ω2∂β = frequencies_from_αβ_derivatives(α, β, model)
    # Chain rule
    ∂Ω1∂a, ∂Ω2∂a = ∂Ω1∂α * ∂α∂a + ∂Ω1∂β * ∂β∂a, ∂Ω2∂α * ∂α∂a + ∂Ω2∂β * ∂β∂a
    ∂Ω1∂e, ∂Ω2∂e = ∂Ω1∂α * ∂α∂e + ∂Ω1∂β * ∂β∂e, ∂Ω2∂α * ∂α∂e + ∂Ω2∂β * ∂β∂e
    return Ω1, Ω2, ∂Ω1∂a, ∂Ω2∂a, ∂Ω1∂e, ∂Ω2∂e
end

########################################################################
#
# (a,e) ↦ frequencies mappings: jacobians
#
########################################################################
"""
    ae_to_αβ_jacobian(a, e, model[, params])

jacobian of the (a,e) ↦ (α,β) mapping, i.e. |∂(α,β)/∂(a,e)|.

@IMPROVE: handle possible of NaN or Inf at this level ?
@IMPROVE: Remove duplicates. Jacobian formula duplicated.
"""
function ae_to_αβ_jacobian(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    # calculate the frequency derivatives
    _, _, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = αβ_from_ae_derivatives(a, e, model, params)
    return abs(∂α∂a*∂β∂e - ∂β∂a*∂α∂e)
end

"""
    ae_to_frequencies_jacobian(a, e, model[, params])

jacobian of the (a,e) ↦ (Ω₁,Ω₂) mapping, i.e. |∂(Ω₁,Ω₂)/∂(a,e)|

@IMPROVE: handle possible of NaN or Inf at this level ?
@IMPROVE: Remove duplicates. Jacobian formula duplicated.
"""
function ae_to_frequencies_jacobian(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    # calculate the frequency derivatives
    _, _, ∂Ω1∂a, ∂Ω2∂a, ∂Ω1∂e, ∂Ω2∂e = frequencies_from_ae_derivatives(a, e, model, params)
    return abs(∂Ω1∂a * ∂Ω2∂e - ∂Ω2∂a * ∂Ω1∂e)
end

########################################################################
#
# frequencies ↦ (a,e) mappings
#
########################################################################
"""
    ae_from_frequencies(Ω1, Ω2, model[, params])
"""
function ae_from_frequencies(
    Ω1::Float64,
    Ω2::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # get the circular orbit (maximum radius) for a given Ω₁,Ω₂. 
    # @IMPROVE: use the stronger constraint.
    # @IMPROVE: use default rmin, rmax (should not matter)
    α = Ω1 / frequency_scale(model)
    acirc = _radius_from_αcircular(α, model, params)
    # and start from ecc=0.5
    # @IMPROVE, is there a more optimal starting eccentricity?
    ainit, einit = acirc, 0.5
    mjac(a::Float64,e::Float64) = frequencies_from_ae_derivatives(a, e, model, params)
    a, e, _, _ = _newton_raphson_ae(ainit, einit, Ω1, Ω2, mjac, params)
    return a, e
end

"""
    ae_from_αβ(α, β, model[, params])
"""
function ae_from_αβ(
    α::Float64,
    β::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # get the circular orbit (maximum radius) for a given α. 
    # @IMPROVE: use the stronger constraint.
    acirc = _radius_from_αcircular(α, model, params)
    # and start from ecc=0.5
    # @IMPROVE, is there a more optimal starting eccentricity?
    ainit, einit = acirc, 0.5
    mjac(a::Float64,e::Float64) = αβ_from_ae_derivatives(a, e, model, params)
    a, e, _, _ = _newton_raphson_ae(ainit, einit, α, β, mjac, params)
    return a, e
end