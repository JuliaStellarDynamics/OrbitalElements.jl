"""

Special treatment for circular orbits where we can use exact relations with potential model
derivatives.
"""


########################################################################
#
# a ↦ circular frequencies
#
########################################################################
"""
    _Ω1circular(r, model[, params])
"""
function _Ω1circular(
    r::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    if r == 0
        return 2 * sqrt(abs(d2ψ(r, model)))
    end

    return sqrt(d2ψ(r, model) + 3 * dψ(r, model) / r)
end

"""
    _Ω2circular(a, model[, params])
"""
function _Ω2circular(
    r::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    if r == 0
        return sqrt(abs(d2ψ(r, model)))
    end

    return sqrt(dψ(r, model) / r)
end

"""
    _αcircular(a, model[, params])

mapping from radius `a` to the dimensionless radial frequency `α` for circular orbits, 
from the epicyclic approximation

`a` stands for the semi-major axis (equivalent to guiding radius r for a circular orbit)
"""
function _αcircular(
    r::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    return _Ω1circular(r, model, params) / frequency_scale(model)
end

"""
    _βcircular(a, model[, params])

mapping from radius `a` to frequency ratio `β` for circular orbits, from the epicyclic
approximation

`a` stands for the semi-major axis (equivalent to guiding radius r for a circular orbit)

Careful treatment implemented for `a == Inf`. (0/0 limit)
@IMPROVE For inner limit, has to be 1/2 for core potentials but not necessarly for cusps.
"""
function _βcircular(
    r::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    if r == 0
        return 1/2
    end
    # At r == Inf, βcircular is 0 / 0 = NaN. (dψ(x)→0, for x→∞)
    # Cure this by estimating the potential derivative growth rate at infinity and 
    # replace by this equivalent.
    # dψ(x) = c x^γ (for x→∞) (⟹ d2ψ(x) = c γ x^(γ-1))
    # ⟹ βcircular(x) = sqrt(c x^γ / x) / sqrt(c γ x^(γ-1) + 3 c x^γ / x)  (for x→∞)
    #                 = 1 / sqrt(3 + γ)
    if r == Inf
        # Estimate growth rate γ of dψ(x)≈x^γ in the outskirts
        rc = radial_scale(model)
        x1, x2 = 1e8 * rc, 1e9 * rc
        γ = round(log(dψ(x2, model) / dψ(x1, model)) / log(x2 / x1))
        if isnan(γ)
            error("Unable to estimate the growth rate of the potential.")
        elseif γ <= -3
            return Inf
        end
        return sqrt(1 / (3 + γ))
    end

    return _Ω2circular(r, model, params) / _Ω1circular(r, model, params)
end


"""
    _β_from_α_circular(α, model[, params])

mapping from dimensionless radial frequency `α to frequency ratio `β` for circular orbits
"""
function _β_from_α_circular(
    α::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    # get the radius corresponding to the circular orbit
    rcirc = _radius_from_αcircular(α, model, params, eps(Float64), eps(Float64))
    return _βcircular(rcirc, model, params)
end


########################################################################
#
# Radius as a function of frequency (mapping inversion)
#
########################################################################
"""
    _radius_from_αcircular(α, model[, params, tolr, tolf])

backwards mapping from `αcircular` to radius

@ASSUMPTION: `αcircular` is a decreasing function of radius
@WARNING: For α larger than the frequency in the center, we do not return an error, while
we probably should.
@IMPROVE: radius is dimensional quantity, we should not be extremizing on it
@IMPROVE: ultimately, tolr and tolf should be inside parameters (maybe inside the 
backward method parameters)
"""
function _radius_from_αcircular(
    α::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters(),
    tolr::Float64=1000.0*eps(Float64),
    tolf::Float64=1000.0*eps(Float64)
)::Float64
    # Edge cases
    if α < _αcircular(Inf, model, params)
        throw(DomainError(α, "Out of bound frequency α"))
        return 0.0
    elseif α > _αcircular(0.0, model, params)
        # @IMPROVE: better be an error (to merge with the other domain error)
        println("Out of bound circular frequency α. Returning r=0 anyway.")
        return 0.0
    elseif α == _αcircular(Inf, model, params)
        return Inf
    elseif α == _αcircular(0.0, model, params)
        return 0.0
    end

    rmin, rmax = 0.0, radial_scale(model)
    # Tweak rmin and rmax to get the frequencies in bounds
    # Safe while loop since we check that the frequency is indeed reachable
    while α > _αcircular(rmin, model, params)
        rmin /= 2
    end
    while α < _αcircular(rmax, model, params)
        rmax *= 2
    end
    # use bisection to find the circular orbit radius corresponding to given frequency
    # @IMPROVE radius is dimensional quantity, we should not be extremizing on it
    return _bisection(
        r -> _αcircular(r, model, params) - α,
        rmin,
        rmax,
        tolx=tolr,
        tolf=tolf
    )
end