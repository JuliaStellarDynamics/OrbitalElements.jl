#####################################
# Abstract types
#####################################
"""
Abstract type for the considered potential
"""
abstract type Potential end

"""
Abstract type for potentials that depend only on r (spherically symmetric or perfectly flat)
"""
abstract type TwoIntegralPotential <: Potential end

"""
Abstract type for potentials that depend on r,z (flattened)
"""
abstract type ThreeIntegralPotential <: Potential end

"""
Abstract type for central potential
"""
abstract type TwoIntegralCentralPotential <: TwoIntegralPotential end
abstract type ThreeIntegralCentralPotential <: ThreeIntegralPotential end

"""
Abstract type for central potential with finite value in the center
"""
abstract type TwoIntegralCentralCorePotential <: TwoIntegralCentralPotential end
abstract type ThreeIntegralCentralCorePotential <: ThreeIntegralCentralPotential end

"""
Abstract type for central potential with infinite value in the center
"""
abstract type TwoIntegralCentralCuspPotential <: TwoIntegralCentralPotential end
abstract type ThreeIntegralCentralCuspPotential <: ThreeIntegralCentralPotential end

#####################################
# Generic functions (not methods)
#####################################
"""
    ψ(r, model::Potential)

Potential value for central potential `model` at radius `r`.
"""
function ψ(r::Number, model::TwoIntegralPotential)
    # ... [model specific implementation] ...
end

"""
    dψ(r, model::Potential)

Potential derivative value for central potential `model` at radius `r`.
"""
function dψ(r::Number, model::TwoIntegralPotential)
    # ... [model specific implementation] ...
end

"""
    d2ψ(r, model::Potential)

Potential second derivative for central potential `model` at radius `r`.
"""
function d2ψ(r::Number, model::TwoIntegralPotential)
    # ... [model specific implementation] ...
end

"""
    frequency_scale(model::Potential)

Frequency scale associated to the potential `model`.
"""
function frequency_scale(model::TwoIntegralPotential)
    # ... [model specific implementation] ...
end

"""
    energy_scale(model::Potential)

Energy scale associated to the potential `model`.
"""
function energy_scale(model::TwoIntegralPotential)
    # ... [model specific implementation] ...
end

"""
    momentum_scale(model::Potential)

Angular momentum scale associated to the potential `model`.
"""
function momentum_scale(model::TwoIntegralPotential)
    # ... [model specific implementation] ...
end

"""
    radial_scale(model::Potential)

Radial scale associated to the potential `model`.
"""
function radial_scale(model::TwoIntegralPotential)
    # ... [model specific implementation] ...
end

#####################################
# Include potentials
#####################################

# two integral potentials
include("isochrone.jl")
include("plummer.jl")
include("hernquist.jl")
include("mestelzang.jl")

# three integral potentials
include("kuzminkutuzov.jl")