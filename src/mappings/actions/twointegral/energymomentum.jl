"""
    (Forward) computations of energy and angular momentum

From semimajor axis and eccentricity to 
Strategies:
-compute energy and angular momentum from the definitions for (a,e)
-compute derivatives of energy and angular momentum
-add switches for near-circular orbits (E,L,dE,dL)
-include Jacobian to transform between (E,L) and (a,e)

"""

########################################################################
#
# (a,e) ↦ (E,L) mapping
#
########################################################################
"""
    EL_from_ae(a, e, model[, params])

energy and angular momentum as a function of semimajor axis and eccentricity
for a given potential `model`.

The generic expressions are defined in ``(r_p,r_a)`` but edges are handled in ``(a,e)``.
See, e.g., equation (65) in Hamilton et al. (2018).

@IMPROVE: right now the edge handling procedure is hard-coded. Should incorparate it in the
OrbitalParameters structure (not only the tolerances)
@IMPROVE: provide stable expressions at infinity
"""
function EL_from_ae(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # Edge cases
    if a == 0 # Values at the center
        return ψ(0., model), 0.
    elseif e == 0 # Values for perfectly circular orbits
        Ecircular = ψ(a, model) + 0.5*a*dψ(a, model)
        Lcircular = sqrt(a)^3 * sqrt(dψ(a, model))
        return  Ecircular, Lcircular
    elseif e == 1 # Values for perfectly radial orbits
        return ψ(2a, model), 0.
    end
    # Handling (a,e)-domain edges through interpolation
    # IMPORTANT : has to be before the generic computation
    fun(atmp::Float64, etmp::Float64) = EL_from_ae(atmp, etmp, model, params)
    res = _interpolate_edges_ae(fun, a, e, params)
    if !isnothing(res)
        return res
    end
    # Generic expression are easily given in (rp,ra)
    rp, ra = rpra_from_ae(a, e)
    E = (ra^2 * ψ(ra, model) - rp^2 * ψ(rp, model)) / (ra^2 - rp^2)
    L = sqrt(2(ψ(ra, model) - ψ(rp, model)) / (rp^(-2) - ra^(-2)))
    return E, L
end

########################################################################
#
# (a,e) ↦ (E,L) mapping: derivatives
#
########################################################################
"""
    EL_from_ae_derivatives(a, e, model[, params])

same as EL_from_ae plus mapping derivatives.

@IMPROVE: right now the edge handling procedure is hard-coded. Should incorparate it in the
OrbitalParameters structure (not only the tolerances)
@IMPROVE: right now the derivative procedure is hard-coded. Should incorparate it 
in the OrbitalParameters structure (not only the steps)
"""
function EL_from_ae_derivatives(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}
    # Function to differentiate
    fun(atemp::Float64, etemp::Float64) = EL_from_ae(atemp, etemp, model, params)
    # Perform differentiation
    floc, ∂f∂a, ∂f∂e = _derivatives_ae(fun, a, e, params)
    # Recast results
    E, L = floc
    ∂E∂a, ∂L∂a = ∂f∂a
    ∂E∂e, ∂L∂e = ∂f∂e

    return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e
end

########################################################################
#
# (a,e) ↦ (E,L) mapping: jacobian
#
########################################################################
"""
    ae_to_EL_jacobian(a, e, model[, params])

jacobian of the (a,e) ↦ (E,L) mapping, i.e. |∂(E,L)/∂(a,e)|.

@IMPROVE: handle possible of NaN or Inf at this level ?
@IMPROVE: Remove duplicates. Jacobian formula duplicated.
"""
function ae_to_EL_jacobian(
    a::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    _, _, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e = EL_from_ae_derivatives(a, e, model, params)
    return abs(∂E∂a * ∂L∂e - ∂L∂a * ∂E∂e)
end

########################################################################
#
# (E,L) ↦ (a,e) mapping
#
########################################################################
"""
    ae_from_EL(E, L, model[, params])
"""
function ae_from_EL(
    E::Float64,
    L::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}
    # Initial guess
    einit = 0.5
    ainitE = try _initial_a_from_E(E, einit, model, params) catch; 0. end
    ainitL = try _initial_a_from_L(L, einit, model, params) catch; 0. end
    ainit = ainitE == 0. && ainitL == 0. ? params.rc : 0.5 * (ainitE + ainitL)
    # Backward mapping
    mjac(a::Float64,e::Float64) = EL_from_ae_derivatives(a, e, model, params)
    a, e, _, _ = _newton_raphson_ae(ainit, einit, E, L, mjac, params)
    return a, e
end

########################################################################
#
# (E,L) ↦ (a,e) mapping: initial guess
#
########################################################################
# Guess using Lcirc poorly behave for eccentricity close to one.
"""
    _initial_a_from_E(E, e, model[, params, tolx, tolf])

perform backwards mapping from energy E for a fixed eccentricity orbit to semi-major axis 
to find a initial guess for semimajor axis.

@ASSUMPTION: important assumption E is a increasing function of semi-major axis 
(at fixed eccentricity)
@IMPROVE: huge duplicates with [`_radius_from_αcircular`](@ref) and other initial 
guess functions
@IMPROVE: ultimately, tolr and tolf should be inside parameters (maybe inside the 
backward method parameters)
"""
function _initial_a_from_E(
    E::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters(),
    tolx::Float64=1000.0*eps(Float64),
    tolf::Float64=1000.0*eps(Float64)
)::Float64
    # check that the input energy is valid
    if E < ψ(0., model)
        error("Too low energy E = ",E)
    elseif E == ψ(0., model)
        return 0.
    elseif E > ψ(Inf, model)
        error("Too high energy E = ",E)
    elseif E == ψ(Inf, model)
        return Inf
    end

    # Tweak xmin and xmax to get the energy in bounds
    # Safe while loop since we check that the energy is indeed reachable
    rmin, rmax = params.rmin, min(params.rmax, 1.e8 * params.rc)
    while ψ(rmin, model) > E
        rmin /= 2
    end
    while ψ(rmax, model) < E
        rmax *= 2
    end

    # use bisection to find the circular orbit radius corresponding to given energy
    # (at fixed eccentricity)
    function rootequation(a::Float64)
        tmpE, _ = EL_from_ae(a, e, model, params)
        return E - tmpE
    end
    return _bisection(rootequation, rmin, rmax, tolx=tolx, tolf=tolf)
end

"""
    _initial_a_from_L(L, e, model[, params, tolx, tolf])

perform backwards mapping from angular momentum for a fixed eccentricity orbit to 
semi-major axis to find a initial guess for semimajor axis.

@IMPROVE: huge duplicates with [`_radius_from_αcircular`](@ref) and other initial 
guess functions
@IMPROVE: ultimately, tolr and tolf should be inside parameters (maybe inside the 
backward method parameters)
"""
function _initial_a_from_L(
    L::Float64,
    e::Float64,
    model::TwoIntegralPotential,
    params::OrbitalParameters=OrbitalParameters(),
    tolx::Float64=1000.0*eps(Float64),
    tolf::Float64=1000.0*eps(Float64)
)::Float64
    # check that the input angular momentum is valid
    if L < 0.
        error("Negative angular momentum L = ",L)
    elseif L == 0.
        return 0.
    elseif L == Inf
        return Inf
    end

    # use bisection to find the circular orbit radius corresponding to given angular
    # momentum (at fixed eccentricity)
    function rootequation(a::Float64)
        _, tmpL = EL_from_ae(a, e, model, params)
        return L - tmpL
    end
    # Tweak xmin and xmax to get the angular momentum in bounds
    # Safe while loop since we check that the angular momentum is indeed reachable
    rmin, rmax = params.rmin, min(params.rmax, 1.e8 * params.rc)
    while rootequation(rmin) < 0 # L(rmin) > objective_L
        rmin /= 2
    end
    while rootequation(rmax) > 0 # L(rmax) < objective_L
        rmax *= 2
    end
    return _bisection(rootequation, rmin, rmax, tolx=tolx, tolf=tolf)
end