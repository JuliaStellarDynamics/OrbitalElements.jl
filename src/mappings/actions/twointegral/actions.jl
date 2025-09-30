

########################################################################
#
# (a,e) ↦ (J,L) mapping
#
########################################################################
"""
    _radial_action_from_ae(a, e, model[, params])

@IMPROVE: right now, Henon anomaly and integrator () are hard-
coded. Should smartly incorparate it in the OrbitalParameters structure.
"""
function _radial_action_from_ae(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    # Edge cases
    if a == 0 || e == 0
        return 0.0
    end
    # Handling (a,e)-domain edges through interpolation
    # IMPORTANT : has to be before the generic computation
    fun(atmp::Float64, etmp::Float64) = _radial_action_from_ae(atmp, etmp, model, params)
    res = _interpolate_edges_ae(fun, a, e, params)
    if !isnothing(res)
        return res
    end
    # Generic computations
    function action_integrand(w::Float64)::Float64
        drdw = radius_from_anomaly_derivative(w, a, e, model, params)
        vrad = radial_velocity(w, a, e, model, params)
        return drdw * vrad
    end
    return (1 / pi) * _integrate_simpson(action_integrand, params.NINT)
end

"""
    actions_from_ae(a, e, model[, params])

mapping from semimajor axis and eccentricity to actions.
"""
function actions_from_ae(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    J = _radial_action_from_ae(a, e, model, params)
    _, L = EL_from_ae(a, e, model, params)
    return J, L
end

########################################################################
#
# (a,e) ↦ (J,L) mapping: derivatives
#
########################################################################
"""
    actions_from_ae_derivatives(a, e, model[, params])

same as actions_from_ae plus mapping derivatives.

@IMPROVE: right now the edge handling procedure is hard-coded. Should incorparate it in the
OrbitalParameters structure (not only the tolerances)
@IMPROVE: right now the derivative procedure is hard-coded. Should incorparate it 
in the OrbitalParameters structure (not only the steps)
"""
function actions_from_ae_derivatives(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}
    # Function to differentiate
    fun(atmp::Float64, etmp::Float64) = actions_from_ae(atmp, etmp, model, params)
    # Perform differentiation
    floc, ∂f∂a, ∂f∂e = _derivatives_ae(fun, a, e, params)
    # Recast results
    J, L = floc
    ∂J∂a, ∂L∂a = ∂f∂a
    ∂J∂e, ∂L∂e = ∂f∂e

    return J, L, ∂J∂a, ∂L∂a, ∂J∂e, ∂L∂e
end

########################################################################
#
# actions ↦ (a,e) mapping
#
########################################################################
"""
    ae_from_actions(J, L, model[, params])
"""
function ae_from_actions(
    J::Float64,
    L::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # Initial guess
    einit = 0.5
    ainitJ = try _initial_a_from_J(J, einit, model, params) catch; 0. end
    ainitL = try _initial_a_from_L(L, einit, model, params) catch; 0. end
    ainit = ainitJ == 0. && ainitL == 0. ? params.rc : 0.5 * (ainitJ + ainitL)
    # Backward mapping
    mjac(a::Float64,e::Float64) = actions_from_ae_derivatives(a, e, model, params)
    a, e, _, _ = _newton_raphson_ae(ainit, einit, J, L, mjac, params)
    return a, e
end

########################################################################
#
# actions ↦ (a,e) mapping: initial guess
#
########################################################################
"""
    _initial_a_from_J(J, e, model[, params, tolx, tolf])

perform backwards mapping from radial action `J` for a fixed eccentricity orbit to 
semi-major axis to find a initial guess for semimajor axis.

@IMPROVE: huge duplicates with [`_radius_from_αcircular`](@ref) and other initial 
guess functions
@IMPROVE: ultimately, tolr and tolf should be inside parameters (maybe inside the 
backward method parameters)
"""
function _initial_a_from_J(
    J::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters(),
    tolx::Float64=1000.0*eps(Float64),
    tolf::Float64=1000.0*eps(Float64)
)::Float64
    # check that the input energy is valid
    if J < 0.
        error("Negative radial action J = ",J)
    elseif J == 0.
        return 0.
    elseif J == Inf
        return Inf
    end

    # use bisection to find the circular orbit radius corresponding to given radial action
    # (at fixed eccentricity)
    rootequation(a::Float64) = J - _radial_action_from_ae(a, e, model, params)
    # Tweak xmin and xmax to get the radial action in bounds
    # Safe while loop since we check that the action is indeed reachable
    rmin, rmax = params.rmin, min(params.rmax, 1.e8 * params.rc)
    while rootequation(rmin) < 0 # J(rmin) > J_objective
        rmin /= 2
    end
    while rootequation(rmax) > 0 # J(rmax) < J_objective
        rmax *= 2
    end

    return _bisection(rootequation, rmin, rmax, tolx=tolx, tolf=tolf)
end
