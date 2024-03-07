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
    _αcircular(a, model)

mapping from radius `a` to the dimensionless radial frequency `α` for circular orbits, 
from the epicyclic approximation

`a` stands for the semi-major axis (equivalent to guiding radius r for a circular orbit)
"""
function _αcircular(a::Float64, model::Potential)::Float64
    if (a == 0.)
        return 2 * sqrt(abs(d2ψ(0., model))) / Ω₀(model)
    end

    return sqrt(d2ψ(a, model) + 3 * dψ(a, model) / a) / Ω₀(model)
end

"""
    _βcircular(a, model)

mapping from radius `a` to frequency ratio `β` for circular orbits, from the epicyclic
approximation

`a` stands for the semi-major axis (equivalent to guiding radius r for a circular orbit)
@IMPROVE For inner limit, has to be 1/2 for core potentials but not necessarly for cusps.
"""
function _βcircular(a::Float64, model::Potential)::Float64
    if (a == 0.)
        return 1/2
    end

    return sqrt(dψ(a, model) / a) / (_αcircular(a, model) * Ω₀(model))
end


"""
    _β_from_α_circular(α, model, params)

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
        min(params.rmax,1.e8*params.rc),
        eps(Float64),
        eps(Float64)
    )

    # @IMPROVE: explain !!
    if rcirc == Inf
        # Estimate growth rate γ of dψ(x)≈x^γ in the outskirts
        x1, x2 = 1e8*params.rc, 1e9*params.rc
        γ = round(log(dψ(x2, model) / dψ(x1, model)) / log(x2 / x1))
        if isnan(γ)
            error("Unable to estimate the growth rate of the potential.")
        elseif γ <= -3
            return Inf
        end
        return sqrt(1 / (3 + γ))
    end

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
        warn("Out of bound circular frequency α. Returning r=0 anyway.")
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