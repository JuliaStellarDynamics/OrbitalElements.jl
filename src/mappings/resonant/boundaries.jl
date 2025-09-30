"""
Boundaries for the resonance variables
"""

########################################################################
#
# (u,v) mapping : ω boundaries (at given resonance number n1, n2)
#
########################################################################
"""
    αminmax(model, params)

maximal and minimal considered radial (dimensionless) frequency

@ASSUMPTION: [`αcircular`](@ref) is a decreasing function of radius
"""
function αminmax(model::TwoIntegralPotential, params::OrbitalParameters)::Tuple{Float64,Float64}
    # Assumption :
    # Ω1circular is a decreasing function of radius
    return _αcircular(params.rmax, model, params), _αcircular(params.rmin, model, params)
end

"""
    frequency_extrema(n1,n2,model,params)

minimal and maximal (dimensionless) "frequency" for a given resonance `ωmin, ωmax = minmax(n⋅Ω/Ω₀)`

@ASSUMPTION: Frequency domain truncated at `αmin` and `αmax` set by `rmax` and `rmin`
@IMPROVE radius is dimensional quantity, we should not be extremizing on it
"""
function frequency_extrema(
    n1::Int64,
    n2::Int64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # For resonance number (0, 0)
    # Avoid infinite loop in `_extremise_noedges` with the flat function `_ωncirc = 0.`
    if n1 == 0 && n2 == 0
        return 0.0, 0.0
    end
    # define the function to extremise, i.e., the adimensional resonant frequency along 
    # the circular line.
    # Works better with n⋅Ω/Ω₀ than n₁α+n₂αβ
    # (+ avoid useless computations in the generic case)
    # @IMPROVE radius is dimensional quantity, we should not be extremizing on it
    function _ωncirc(r::Float64)::Float64
        return (
            (n1 * _Ω1circular(r, model, params) + n2 * _Ω2circular(r, model, params))
            / frequency_scale(model)
        )
    end
    # If rmax is infinite, bisection search on a bounded interval
    rc, rmin, rmax = params.rc, params.rmin, params.rmax
    xext = _extremise_noedges(_ωncirc, rmin, min(rmax, 1e8 * rc))
    # The extreme values of n.Ω/Ω₀ is either :
    #   - on the radial line, at α = αmin or αmax
    #   - along the circular velocity (extreme α included)
    αmin, αmax = αminmax(model, params)
    ωmin = min(
        _ωncirc(xext),
        _ωncirc(rmin),
        _ωncirc(rmax), 
        (n1 + n2 / 2) * αmin,
        (n1 + n2 / 2) * αmax
    )
    ωmax = max(
        _ωncirc(xext),
        _ωncirc(rmin),
        _ωncirc(rmax),
        (n1 + n2 / 2) * αmin,
        (n1 + n2 / 2) * αmax
    )

    return ωmin, ωmax
end

########################################################################
#
# Resonance struct to store the resonance number and frequencies extrema
#
########################################################################
struct Resonance
    number::Tuple{Int64,Int64}
    frequency_extrema::Tuple{Float64,Float64}
end
"""
    Resonance(n1, n2, ωmin, ωmax)

Structure to store the resonance number and its associated extremal frequencies.

# Fields
- `number::Tuple{Int64,Int64}`
- `frequency_extrema::Tuple{Float64,Float64}`
"""
function Resonance(n1::Int64, n2::Int64, ωmin::Float64, ωmax::Float64)
    return Resonance((n1, n2), (ωmin, ωmax))
end
function Resonance(
    n1::Int64, 
    n2::Int64, 
    model::TwoIntegralPotential, 
    params::OrbitalParameters=OrbitalParameters()
)
    ωmin, ωmax = frequency_extrema(n1, n2, model, params)
    return Resonance((n1, n2), (ωmin, ωmax))
end

frequency_extrema(res::Resonance) = res.frequency_extrema
resonance_number(res::Resonance) = res.number

########################################################################
#
# (u,v) mapping : v boundary (at given u, n1, n2)
#
########################################################################
"""
    v_boundaries(u, res, model, params)

for a given resonance, at a specific value of u, find the v coordinate boundaries.
"""
function v_boundaries(
    u::Float64,
    res::Resonance,
    model::TwoIntegralCentralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    rmin, rmax = params.rmin, params.rmax
    αmin, αmax = αminmax(model, params)

    # ωn(u) : value of the resonance line
    ωmin, ωmax = frequency_extrema(res)
    hval = _ωn(u, ωmin, ωmax)

    # βcircular as a function of αcircular
    βc(αc::Float64)::Float64 = _β_from_α_circular(αc, model, params)

    n1, n2 = resonance_number(res)
    if n2 == 0 # v = β
        #####
        # (B9) Fouvry & Prunet : 1st inequality
        #####
        vmin = 0.5
        #####
        # (B9) Fouvry & Prunet : 2nd inequality
        #####
        vmax = n1 == 0 ? vmin : βc(hval/n1)
    else # v = α

        #####
        # (B10) Fouvry & Prunet : 1st & 2nd inequalities
        #####
        vmin = αmin
        vmax = αmax

        #####
        # Constraint for hval = 0.
        # Line of constant β (or v = 0.)
        #####
        # β = - n1/n2 ≥ 1/2
        if hval == 0 && -n1/n2 < 1/2
            vmin = 0.
            vmax = 0.
        end
        # β = - n1/n2 ≤ βc(v)
        # Hypothesis : βc is a decreasing function of α (=v)
        if hval == 0 && -n1/n2 >= βc(αmin)
            vmin = αmin
            vmax = αmin
        end

        #####
        # (B10) Fouvry & Prunet : 3rd inequality
        #####
        # Hiting radial boundary β = 1/2
        # hval = (n1+n2/2)v
        radon = n1 + n2 / 2 # Radial orbit equivalent n
        if n2 * hval > 0 && radon * hval > 0
            vmax = min(vmax, hval/radon) # Updating vmax
        elseif n2 * hval < 0 && radon * hval > 0
            vmin = max(vmin, hval/radon) # Updating vmin
        end

        #####
        # (B10) Fouvry & Prunet : 4th inequality
        #####
        # boolean monotonic stands for x -> ωncirc(x) monotonicity
        if n1 * n2 > 0 # Do not search for vbound
            monotonic = true
        else
            # First look for vbound in the asked boundary
            vbound = _α_inner_extremum(res, rmin, rmax, model, params)

            # Extreme boundary to look for vbound
            locrmin, locrmax = 0.0, Inf

            if (vbound != αmin) && (vbound != αmax)
                monotonic = false
            elseif (rmin > locrmin) || (rmax < locrmax)
                # If vbound not in the asked boundary
                # verify that it should indeed not exist
                locrmin, locrmax = min(rmin, locrmin), max(rmax, locrmax)
                vbound = _α_inner_extremum(res, locrmin, locrmax, model, params)
                # If we have found an extremum in the inner domain
                if !(vbound in (
                    _αcircular(locrmin, model, params), 
                    _αcircular(locrmax, model, params)
                    )
                )
                    monotonic = false
                else 
                    monotonic = true
                end
            else
                monotonic = true
            end
        end

        # Constraint equation
        rootequation(v::Float64)::Float64 = hval - n1 * v - n2 * v * βc(v)

        if monotonic # ωncirc(x) is monotonic
            vlim = try _bisection(rootequation, αmin, αmax) catch; -1 end
            if vlim != -1 && hval * n2 > 0
                vmin = max(vmin, vlim)
            elseif vlim != -1 && hval * n2 < 0
                vmax = min(vmax, vlim)
            end
        else # ωncirc(x) is not monotonic
            # Search crossing in [αmin,vbound]
            vmin2 = try _bisection(rootequation, αmin, vbound) catch; -1 end
            if (vmin2 != -1)
                vmin = max(vmin, vmin2)
            end
            # Search crossing in [vbound,αmax]
            vmax2 = try _bisection(rootequation, vbound, αmax) catch; -1 end
            if (vmax2 != -1)
                vmax = min(vmax, vmax2)
            end
        end
    end

    return vmin, vmax
end


"""
    _ωn(u,ωmin,ωmax)

dimensionless frequency of the resonance line `u`.

return h_n(u) = ω_n(u), a helper quantity
See, equation (B8) in Fouvry & Prunet (2022)
"""
function _ωn(u::Float64, ωmin::Float64, ωmax::Float64)::Float64
    if u == 1
        return ωmax
    elseif u == -1
        return ωmin
    end 
    return (ωmax + ωmin + u * (ωmax - ωmin)) / 2
end


"""
    _α_inner_extremum(res, rmin, rmax, model[, params])

location `α` of the extremum of the resonance frequency along the circular line limited to 
`]rmin,rmax[`.

Importantly allows to search for rmin and rmax different from the ones in params.

@IMPROVE radius is dimensional quantity, we should not be extremizing on it
"""
function _α_inner_extremum(
    res::Resonance,
    rmin::Float64,
    rmax::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    # define the function to extremise, i.e., the adimensional resonant frequency along 
    # the circular line.
    # Works better with n⋅Ω/Ω₀ than n₁α+n₂αβ
    # (+ avoid useless computations in the generic case)
    n1, n2 = resonance_number(res)
    # @IMPROVE radius is dimensional quantity, we should not be extremizing on it
    function _ωncirc(r::Float64)::Float64
        return (
            (n1 * _Ω1circular(r, model, params) + n2 * _Ω2circular(r, model, params))
            / frequency_scale(model)
        )
    end
    # If rmax is infinite, bisection search on a bounded interval
    locrmax = min(rmax, 1e8 * radial_scale(model))
    xext = _extremise_noedges(_ωncirc, rmin, locrmax)

    # If the extremum is reached at the imposed maximal boundary
    # Use the true rmax (not the artificial 1e8 * rc, which is here to handle Inf)
    if xext == locrmax
       xext = rmax
    end

    return _αcircular(xext, model, params)
end
