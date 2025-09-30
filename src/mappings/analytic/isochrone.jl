


########################################################################
#
# (a,e) ↔ (E,L) mappings
#
########################################################################
"""
for isochrone analytical version, see equation (G8) in Fouvry+21
"""
function EL_from_ae(
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    rp, ra = rpra_from_ae(a, e)
    xp, xa = (rp, ra) ./ radial_scale(model) # dimensionless radius
    sp, sa = _spsa_from_rpra(rp, ra, model, params) # extremal anomalies on the orbit
    Ẽ = 1 / (sp + sa) # dimensionless energy
    L̃ = sqrt(2) * xp * xa / sqrt((1 + sp) * (1 + sa) * (sp + sa)) # dimensionless momentum
    return energy_scale(model) * Ẽ, momentum_scale(model) * L̃ 
end

"""    
for isochrone analytical version, 
see equations @TOCOMPLETE in Fouvry+21
@IMPROVE, uses a floor to avoid any sqrt problems with circular orbits
@IMPROVE, add domain check and border treatment (`E==0``)
"""
function ae_from_EL(
    E::Float64,
    L::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    Ẽ = E / energy_scale(model) # dimensionless energy
    xc = 1 / (2Ẽ) - 1
    eccov = sqrt(max(0.0, 1 - (L / momentum_scale(model))^2 * (1 / xc) * (1 + (1 / xc))))
    xp = sqrt(max(0.0,(2 + xc * (1 - eccov)) * (xc * (1 - eccov))))
    xa = sqrt(max(0.0,(2 + xc * (1 + eccov)) * (xc * (1 + eccov))))
    rp, ra = radial_scale(model) .* (xp, xa) # dimensional peri and apocentres
    return ae_from_rpra(rp, ra)
end

"""
for isochrone analytical version
@IMPROVE: Define it for the isochrone !
"""
# function EL_from_ae_derivatives(
#     a::Float64,
#     e::Float64,
#     model::AnalyticIsochrone,
#     params::OrbitalParameters=OrbitalParameters()
# )
#     error("not analytically defined for isochrone")
# end

########################################################################
#
# (a,e) ↦ (J,L) mappings
#
########################################################################
"""
for isochrone analytical version, see equation (G3) in Fouvry+21
@IMPROVE: use more explicit dimensionalisation.
"""
function _radial_action_from_ae(
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    E, L = EL_from_ae(a, e, model, params)
    G, M, bc = model.G, model.M, model.bc
    return G * M / sqrt(-2E) - (L + sqrt(L^2 + 4G * M * bc)) / 2
end

########################################################################
#
# (a,e) ↔ frequencies mappings: integrands
#
########################################################################
"""
the analytic expression for Theta for the isochrone profile are given in 
Fouvry, Hamilton, Rozier, Pichon eq. (G10). It has the major advantage of 
always being well-posed for -1<u<1
"""
function _Θ(
    w::Float64,
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    r = radius_from_anomaly(w, a, e, model, params)
    rp, ra = rpra_from_ae(a, e)
    # dimensionless pericentre, apocentre, and radius
    x, xp, xa = (r, rp, ra) ./ radial_scale(model)
    s = _s_from_r(r, model, params)
    sp, sa = _spsa_from_rpra(rp, ra, model, params)

    # Analytical expression of (dr/du)(1/vr), that is always well-posed
    return (
        3x * sqrt(
            (s + sp) 
            * (s + sa) 
            * (sp + sa)
            / (x + xp) 
            / (x + xa) 
            / (4 - w^2)
        )
        / (sqrt(2) * frequency_scale(model))
    )
end

########################################################################
#
# (a, e) ↔ frequencies mappings: circular
#
########################################################################
function _Ω1circular(
    r::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    s = _s_from_r(r, model, params)
    return frequency_scale(model) / s^(3/2)
end

function _Ω2circular(
    r::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    s = _s_from_r(r, model, params)
    return frequency_scale(model) / (sqrt(s) * (1 + s))
end

function _βcircular(
    r::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    s = _s_from_r(r, model, params)
    # Cure at a == Inf, βcircular is 0 / 0 = NaN.
    if r == Inf
        return one(a)
    end

    return s / (1 + s)
end

"""
for isochrone analytical version, see equation (C6) in Fouvry & Prunet (2022)
"""
function _β_from_α_circular(
    α::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    # Circular orbits: s = α^(-2/3)
    # and 1 / (1 + 1/x) = x / (1 + x)
    return 1/(1 + α^(2/3))
end

function _radius_from_αcircular(
    α::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters(),
    tolr::Float64=1000.0*eps(Float64),
    tolf::Float64=1000.0*eps(Float64)
)::Float64
    # Circular orbits: α = s^(-3/2)
    s = α^(-2/3)
    return _r_from_s(s, model, params)
end

########################################################################
#
# (a, e) ↔ frequencies mappings
#
########################################################################
"""
for isochrone analytical version, see equations (G5-G7) in Fouvry+21
"""
function αβ_from_ae(
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    rp, ra = rpra_from_ae(a, e)
    xp, xa = (rp, ra) ./ radial_scale(model)
    sp, sa = _spsa_from_rpra(rp, ra, model, params)
    return (2 / (sp + sa))^(3/2), (1 + (xp * xa) / ((1 + sp) * (1 + sa))) / 2
end

function ae_from_αβ(
    α::Float64,
    β::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    # (α,β) ↦ (E,L), analytic expressions specific to isochrone potential
    E, L = _EL_from_αβ(α, β, model, params)
    # (E,L) ↦ (a,e), analytic expressions specific to isochrone potential
    a, e = ae_from_EL(E, L, model, params)
    return a, e
end

# `frequencies_from_ae(a, e, ::AnalyticIsochrone)`
# does not need to be (re)defined. 
# automatically uses `frequencies_from_αβ(αβ_from_ae(a, e, ::AnalyticIsochrone),...)`

# however, `ae_from_frequencies` does, because the generic one does not use 
# `αβ_from_frequencies` but rather try its own inversion (with its own initialization).
function ae_from_frequencies(
    Ω1::Float64,
    Ω2::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    α, β = αβ_from_frequencies(Ω1, Ω2, model)
    return ae_from_αβ(α, β, model, params)
end

function ae_to_αβ_jacobian(
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    # @IMPROVE: uses the generic (non-analytic) version for now
    jacaetoEL = ae_to_EL_jacobian(a, e, model, params)
    E, L = EL_from_ae(a, e, model, params)
    jacELtoαβ = _EL_to_αβ_jacobian(E, L, model, params)
    # Chain rule
    return jacaetoEL * jacELtoαβ
end

function ae_to_frequencies_jacobian(
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    jacaetoαβ = ae_to_αβ_jacobian(a, e, model, params)
    α, β = αβ_from_ae(a, e, model, params)
    jacαβtofrequencies = αβ_to_frequencies_jacobian(α, β, model)
    # Chain rule
    return jacaetoαβ * jacαβtofrequencies
end


########################################################################
#
# (E,L) ↔ frequencies mappings
# @WARNING: This is specific to isochrone, not defined for other (non-analytic) potentials. 
# However we use them inside the (exported) analytic mappings
# We do not export these functions! 
#
########################################################################
"""
analytical isochrone energy and angular momentum as a function of frequency ratios
see equations @TOCOMPLETE in Fouvry+21
"""
function _EL_from_αβ(
    α::Float64,
    β::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    return (
        energy_scale(model) * α^(2/3) / 2, 
        momentum_scale(model) * (2β - 1) / (sqrt(β * (1 - β)))
    )
end

"""
analytical isochrone frequency ratios as a function of energy and angular momentum
"""
function _αβ_from_EL(
    E::Float64,
    L::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    # Dimensionless energy and momentum
    Ẽ, L̃ = E / energy_scale(model), L / momentum_scale(model) 
    return 2Ẽ^(3/2), (1 + 1 / sqrt(1 + 4 / L̃^2)) / 2
end

"""
analytic jacobian of the (α,β) ↦ (E,L) mapping, i.e. |∂(E,L)/∂(α,β)|.
"""
function _αβ_to_EL_jacobian(
    α::Float64,
    β::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    return abs(
        energy_scale(model) * momentum_scale(model) 
        / (6 * α^(1/3) * (β * (1 - β))^(3/2))
    )
end

"""
analytic jacobian of the (E,L) ↦ (α,β) mapping, i.e. |∂(α,β)/∂(E,L)|.
"""
function _EL_to_αβ_jacobian(
    E::Float64,
    L::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    Ẽ, L̃ = E / energy_scale(model), L / momentum_scale(model)
    return abs(
        6 * sqrt(2Ẽ)
        / (
            energy_scale(model) 
            * momentum_scale(model)
            * (1 + 4 / L̃^2)^(3/2)
            * L̃^3
        )
    )
end

########################################################################
#
# (α, β) ↔ (u,v) mappings
#
########################################################################
"""
for isochrone analytical version, only the cut in `αmin` (i.e., `params.rmax`) is used
@IMPROVE add the cut in the centre (`αmax`, i.e., `params.rmin`)
"""
function frequency_extrema(
    n1::Int64,
    n2::Int64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    αmin, _ = αminmax(model, params) # Minimum value allowed for α
    omg1 = (n1 + n2 / 2) * αmin # Bottom left value of omega
    omg2 = (n1 + n2 * _β_from_α_circular(αmin, model)) * αmin # Top left value of omega
    omg3 = (n1 + n2 / 2) * _αcircular(0., model) # Bottom right value of omega
    #####
    omgmin, omgmax = extrema((omg1, omg2, omg3)) # Computing the min/max values from the edge
    #####
    if n1 == 0 || n2 == 0 # In these cases, min/max are reached in the edges
        return omgmin, omgmax # Output
    else # In other cases, there might be a maximum along the circular orbits
        Delta = n2 * (n2 - 24 * n1) # Value of the discriminant of P(q)
        if Delta <= 0 # There are no maximum along circular orbits
            return omgmin, omgmax
        else # The discriminant is positive, there might be an extremum along circular orbits
            qcut = 1.0 + αmin^(2/3) # Minimum value of q
            sqDelta = sqrt(Delta) # Square root of the discriminant
            qm = (-n2 - sqDelta)/(6 * n1) # First  root of P(q)
            qp = (-n2 + sqDelta)/(6 * n1) # Second root of P(q)
            #####
            # We sort the roots. At this stage, we necessarily have qmin <= 0.0
            _, qmax = minmax(qm,qp)
            #####
            # We have reached an extremum along circular orbits in the allowed domain
            if qcut <= qmax <= 2
                # Value of α for which an extremum is reached along circular orbits
                αcirc = (qmax - 1)^(3/2)
                # Value of the extremum
                omgcirc = (n1 + n2 * _β_from_α_circular(αcirc, model)) * αcirc
                #####
                # Determining the new minimum and maximum values of the resonance frequency
                omgmin, omgmax = extrema((omgmin, omgmax, omgcirc)) 
                return omgmin, omgmax # Output
            else # We do not reach an extremum along circular orbits in the allowed domain
                return omgmin, omgmax # Output
            end
        end
    end
end

"""
for isochrone analytical version, only the cut in `αmin` (i.e., `params.rmax`) is used
@IMPROVE add the cut in the centre (`αmax`, i.e., `params.rmin`). Caution with extrema 
search: probably not just replacing `αcentre` by `αmax`
@IMPROVE: Lot of code duplicates in this function.
@IMPROVE: Does not work on some edge cases (e.g. res=(-5,-5) and u=-1)
"""
function v_boundaries(
    u::Float64,
    res::Resonance,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # Constants along the function
    n1, n2 = res.number
    αmin, _ = αminmax(model, params) # Minimum value allowed for α
    # Maximal value allowed for α
    # @IMPROVE: use αmax given by `αminmax` with `params.rmin`
    αcentre = _αcircular(0.0, model)
    omgmin, omgmax = res.frequency_extrema
    exph = 0.5*(omgmax + omgmin + u*(omgmax - omgmin)) # Expression of h(u)
    # Resonant frequency along the circular line
    function _ωn_circular(α::Float64)
        return (n1 + n2 * _β_from_α_circular(α, model, params)) * α 
    end
    #####
    # For n2=0, the boundary takes a simple form
    if n2 == 0
        vm = 0.5 # Minimum bound
        vp = n1 == 0 ? vm : _β_from_α_circular(exph/n1, model) # Maximum bound
        return vm, vp # Output
    end
    #####
    # For n2!=0, the boundary is more complicated
    #####
    # First, we account for the easy constraints
    vm = αmin # Minimum bound
    vp = αcentre # Maximum bound
    #####
    kappa = n1 + n2 / 2 # Definition of the constant kappa
    ####
    # Accounting for the third constraint
    ####
    if n2 < 0
        if kappa < 0
            vp = min(vp, exph/kappa) # Updating vp
        elseif kappa > 0
            vm = max(vm, exph/kappa) # Updating vm
        end
    elseif n2 > 0
        if kappa < 0
            vm = max(vm, exph/kappa) # Updating vm
        elseif kappa > 0
            vp = min(vp, exph/kappa) # Updating vp
        end
    end
    #####
    # Accounting for the fourth constraint
    #####
    # Values on the left/right of the interval
    omgl = _ωn_circular(αmin) - exph
    omgr = _ωn_circular(αcentre) - exph
    # The case n1=0 is easy to account for
    if n1 == 0
        if omgl * omgr < 0.0 # We can search for a root
            # Finding the transition coordinate
            vlim = _bisection(v -> _ωn_circular(v) - exph, αmin, αcentre)
            vm = max(vm, vlim) # Updating vm
        end
        #####
        return vm, vp # Output
    end
    #####
    # Now, we know for sure that n1!=0 and n2!=0
    Delta = n2 * (n2 - 24 * n1) # Value of the discriminant
    #####
    if Delta <= 0 # There is no extrema along circular orbits
        if omgl * omgr < 0.0 # We can search for a root
            # Finding the transition coordinate
            vlim = _bisection(v -> _ωn_circular(v) - exph, αmin, αcentre)
            if n1 * n2 > 0
                vm = max(vm, vlim) # Updating vm
            elseif n1 * n2 < 0
                vp = min(vp, vlim) # Updating vp
            end
        end

        return vm, vp # Output
    else # Otherwise, we have Delta > 0, and we might find an extremum along circular orbits
        qc = 1.0 + αmin^(2/3) # Left value of the range in q
        #####
        qm = (-n2 - sqrt(Delta)) / (6 * n1) # Value of q{-}
        qp = (-n2 + sqrt(Delta)) / (6 * n1) # Value of q{+}
        #####
        _, qmax = minmax(qm, qp) # Sorting the two roots
        #####
        # The considered range q \in [qc,2] is on the right of the two roots of P(q)
        if qmax <= qc
            if omgl * omgr < 0.0 # We can search for a root
                # Finding the transition coordinate
                vlim = _bisection(v -> _ωn_circular(v) - exph, αmin, αcentre)
                if n1*n2 > 0
                    vm = max(vm, vlim) # Updating vm
                elseif n1*n2 < 0
                    vp = min(vp, vlim) # Updating vp
                end
            end

            return vm, vp # Output
        #####
        # The considered range q \in [qc,2] is within the two roots of P(q)
        # The only change with the above case is a change in the sign
        # for the condition on n1*n2
        elseif 2.0 <= qmax
            if omgl * omgr < 0.0 # We can search for a root
                # Finding the transition coordinate
                vlim = _bisection(v -> _ωn_circular(v) - exph, αmin, αcentre)
                if n1 * n2 > 0
                    vp = min(vp, vlim) # Updating vp
                elseif n1 * n2 < 0
                    vm = max(vm, vlim) # Updating vm
                end
            end
            return vm, vp
        #####
        # The considered range q \in [qc,2] contains qmax,
        # i.e. contains a root of P(q)
        # This is the case that is more difficult to consider
        else
            vmax = (qmax - 1.0)^(3/2) # Position in v where the monotonicity changes
            #####
            # Value of _ωn_circular at critical points of the interval 
            # (boundaries and extrema)
            omgLEFT = _ωn_circular(αmin)
            omgMIDD = _ωn_circular(vmax)
            omgRGHT = _ωn_circular(αcentre)
            #####
            if n1 < 0 && n2 < 0
                #####
                error("BUG -- getvmvp: CASE 1.")
                #####
            elseif n1 < 0 && n2 > 0
                #####
                if omgLEFT <= exph <= omgMIDD
                    # Finding the transition
                    vlim = _bisection(v -> _ωn_circular(v) - exph, αmin, vmax)
                    vm = max(vm, vlim) # Updating vm
                end
                if omgRGHT <= exph <= omgMIDD
                    # Finding the transition
                    vlim = _bisection(v -> _ωn_circular(v) - exph, vmax, αcentre)
                    vp = min(vp,vlim) # Updating vp
                end
            elseif n1 > 0 && n2 < 0
                if omgMIDD <= exph <= omgLEFT
                    # Finding the transition
                    vlim = _bisection(v -> _ωn_circular(v) - exph, αmin, vmax)
                    vm = max(vm, vlim) # Updating vm
                end
                if omgMIDD <= exph <= omgRGHT
                    # Finding the transition
                    vlim = _bisection(v -> _ωn_circular(v) - exph, vmax, αcentre)
                    vp = min(vp, vlim) # Updating vp
                end
            elseif n1 > 0 && n2 > 0
                error("BUG -- getvmvp: Case 2.")
            end

            return vm, vp
        end
    end
end