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
    _Ω1circular(a, model)
"""
function _Ω1circular(a::Float64, model::Potential)::Float64
    if a == 0.
        return 2 * sqrt(abs(d2ψ(0., model)))
    end

    return sqrt(d2ψ(a, model) + 3 * dψ(a, model) / a)
end

"""
    _Ω2circular(a, model)
"""
function _Ω2circular(a::Float64, model::Potential)::Float64
    if (a == 0.)
        return sqrt(abs(d2ψ(0., model)))
    end

    return sqrt(dψ(a, model)/a)
end

"""
    _αcircular(a, model)

mapping from radius `a` to the dimensionless radial frequency `α` for circular orbits, 
from the epicyclic approximation

`a` stands for the semi-major axis (equivalent to guiding radius r for a circular orbit)
"""
function _αcircular(a::Float64, model::Potential)::Float64
    return _Ω1circular(a, model) / frequency_scale(model)
end

"""
    _βcircular(a, model)

mapping from radius `a` to frequency ratio `β` for circular orbits, from the epicyclic
approximation

`a` stands for the semi-major axis (equivalent to guiding radius r for a circular orbit)

Careful treatment implemented for `a == Inf`. (0/0 limit)
@IMPROVE For inner limit, has to be 1/2 for core potentials but not necessarly for cusps.
"""
function _βcircular(a::Float64, model::Potential)::Float64
    if a == 0.
        return 1/2
    end
    # At a == Inf, βcircular is 0 / 0 = NaN. (dψ(x)→0, for x→∞)
    # Cure this by estimating the potential derivative growth rate at infinity and 
    # replace by this equivalent.
    # dψ(x) = c x^γ (for x→∞) (⟹ d2ψ(x) = c γ x^(γ-1))
    # ⟹ βcircular(x) = sqrt(c x^γ / x) / sqrt(c γ x^(γ-1) + 3 c x^γ / x)  (for x→∞)
    #                 = 
    if a == Inf
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

    return _Ω2circular(a, model) / _Ω1circular(a, model)
end


"""
    _β_from_α_circular(α, model[, params])

mapping from dimensionless radial frequency `α to frequency ratio `β` for circular orbits
"""
function _β_from_α_circular(
    α::Float64,
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    # get the radius corresponding to the circular orbit
    rcirc = _radius_from_αcircular(
        α,
        model,
        params.rmin,
        min(params.rmax, 1.e8 * params.rc),
        eps(Float64),
        eps(Float64)
    )

    return _βcircular(rcirc, model)
end


########################################################################
#
# Radius as a function of frequency (mapping inversion)
#
########################################################################
"""
    _radius_from_αcircular(α, model[, rmin, rmax, tolx, tolf])

backwards mapping from `αcircular` to radius

can tune [rmin,rmax] for extra optimisation (but not needed)
@ASSUMPTION: `αcircular` is a decreasing function of radius
@WARNING: For α larger than the frequency in the center, we do not return an error, while
we probably should.
"""
function _radius_from_αcircular(
    α::Float64,
    model::Potential,
    rmin::Float64=0.,
    rmax::Float64=1.,
    tolr::Float64=1000.0*eps(Float64),
    tolf::Float64=1000.0*eps(Float64)
)::Float64
    # Edge cases
    if α < _αcircular(Inf, model)
        throw(DomainError(α, "Out of bound frequency α"))
        return 0.
    elseif α > _αcircular(0., model)
        # @IMPROVE: better be an error (to merge with the other domain error)
        println("Out of bound circular frequency α. Returning r=0 anyway.")
        return 0.
    elseif α == _αcircular(Inf, model)
        return Inf
    elseif α == _αcircular(0., model)
        return 0.
    end

    # Tweak xmin and xmax to get the frequencies in bounds
    # Safe while loop since we check that the frequency is indeed reachable
    while α > _αcircular(rmin, model)
        rmin /= 2
    end
    while α < _αcircular(rmax, model)
        rmax *= 2
    end
    # use bisection to find the circular orbit radius corresponding to given frequency
    return _bisection(r -> _αcircular(r, model) - α, rmin, rmax, tolx=tolr, tolf=tolf)
end