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
    model::Potential
)::Tuple{Float64,Float64}
    return αβ_from_frequencies(Ω1, Ω2, Ω₀(model))
end


"""
    frequencies_from_αβ(α, β, model)

frequencies from frequency ratios
"""
function frequencies_from_αβ(α::Float64, β::Float64, Ω0::Float64)::Tuple{Float64,Float64}
    return Ω0*α, Ω0*α*β
end
function frequencies_from_αβ(
    Ω1::Float64,
    Ω2::Float64,
    model::Potential
)::Tuple{Float64,Float64}
    return frequencies_from_αβ(Ω1, Ω2, Ω₀(model))
end

"""
    _dfrequencies_from_dαβ(α, β, dα, dβ, model)

convert frequencies ratios derivatives to frequencies derivatives.
"""
function _dfrequencies_from_dαβ(
    α::Float64,
    β::Float64,
    dα::Float64,
    dβ::Float64,
    Ω0::Float64
)::Tuple{Float64,Float64}
    return Ω0 * dα, Ω0 * (dα * β + α * dβ)
end
function _dfrequencies_from_dαβ(
    α::Float64,
    β::Float64,
    dα::Float64,
    dβ::Float64,
    model::Potential
)::Tuple{Float64,Float64}
    return _dfrequencies_from_dαβ(α, β, dα, dβ, Ω₀(model))
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
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # Edge cases
    if a == 0 || e == 0 # Inner or circular frequencies
        return _αcircular(a, model), _βcircular(a, model)
    end
    # Handling (a,e)-domain edges through interpolation
    # IMPORTANT : has to be before the generic computation
    fun(atmp::Float64, etmp::Float64) = αβ_from_ae(atmp, etmp, model, params)
    res = _interpolate_edges_ae(fun, a, e, params)
    if !isnothing(res)
        return res
    end
    # Generic computations
    # 1/α = Ω₀ / π * \int_{-1}^1 du Θ(u); β = L / π * \int_{-1}^1 du Θ(u) / [r(u)]^2
    # using Θ calculations to compute frequencies: leans heavily on Θ from 
    # frequencies/anomaly.jl
    # @IMPROVE: EDGE could be adaptive based on circularity and small-ness of rperi
    # currently using Simpson's 1/3 rule only (hard-coded)
    # @IMPROVE: free the integration scheme, i.e., make it a parameter
    function invαβ_integrands(u::Float64)::Tuple{Float64,Float64}
        # push integration forward on two different quantities: Θ(u),Θ(u)/r^2(u)
        integrand = Θ(u,a,e,model,params)
        return integrand, integrand / (radius_from_anomaly(u, a, e)^2)
    end
    accum1, accum2 = _integrate_simpson(invαβ_integrands, params.NINT)
    _, L = EL_from_ae(a, e, model, params)
    α = pi / (Ω₀(model) * accum1)
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
    model::Potential,
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

same as αβ_from_ae plus mapping derivatives.

@IMPROVE: right now the edge handling procedure is hard-coded. Should incorparate it in the
OrbitalParameters structure (not only the tolerances)
@IMPROVE: right now the derivative procedure is hard-coded. Should incorparate it 
in the OrbitalParameters structure (not only the steps)
"""
function αβ_from_ae_derivatives(
    a::Float64,
    e::Float64,
    model::Potential,
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
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}
    α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = αβ_from_ae_derivatives(a, e, model, params)

    Ω1, Ω2 = frequencies_from_αβ(α, β, model)
    ∂Ω1∂a, ∂Ω2∂a = _dfrequencies_from_dαβ(α, β, ∂α∂a, ∂β∂a, model)
    ∂Ω1∂e, ∂Ω2∂e = _dfrequencies_from_dαβ(α, β, ∂α∂e, ∂β∂e, model)

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
    model::Potential,
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
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    # calculate the frequency derivatives
    _, _, ∂Ω1∂a, ∂Ω2∂a, ∂Ω1∂e, ∂Ω2∂e = frequencies_from_ae_derivatives(model, a, e, params)
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
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # get the circular orbit (maximum radius) for a given Ω₁,Ω₂. 
    # @IMPROVE: use the stronger constraint.
    # @IMPROVE: use default rmin, rmax (should not matter)
    acirc = _radius_from_αcircular(
        Ω1/Ω₀(model), 
        model, 
        params.rmin, 
        min(params.rmax, 1e8*params.rc)
    )
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
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # get the circular orbit (maximum radius) for a given Ω₁,Ω₂. 
    # @IMPROVE: use the stronger constraint.
    # @IMPROVE: use default rmin, rmax (should not matter)
    acirc = _radius_from_αcircular(α, model, params.rmin, min(params.rmax, 1e8*params.rc))
    # and start from ecc=0.5
    # @IMPROVE, is there a more optimal starting eccentricity?
    ainit, einit = acirc, 0.5
    mjac(a::Float64,e::Float64) = αβ_from_ae_derivatives(a, e, model, params)
    a, e, _, _ = _newton_raphson_ae(ainit, einit, α, β, mjac, params)
    return a, e
end