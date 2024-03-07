#####################################
# Abstract types
#####################################
"""
Abstract type for the considered potential
"""
abstract type Potential end

"""
Abstract type for central potential
"""
abstract type CentralPotential <: Potential end

"""
Abstract type for central potential with finite value in the center
"""
abstract type CentralCorePotential <: CentralPotential end

"""
Abstract type for central potential with infinite value in the center
"""
abstract type CentralCuspPotential <: CentralPotential end

#####################################
# Generic functions (not methods)
#####################################
"""
    ψ(r, model::Potential)

Potential value for central potential `model` at radius `r`.
"""
function ψ(r::Number, model::Potential)
    # ... [model specific implementation] ...
end

"""
    dψ(r, model::Potential)

Potential derivative value for central potential `model` at radius `r`.
"""
function dψ(r::Number, model::Potential)
    # ... [model specific implementation] ...
end

"""
    d2ψ(r, model::Potential)

Potential second derivative for central potential `model` at radius `r`.
"""
function d2ψ(r::Number, model::Potential)
    # ... [model specific implementation] ...
end

"""
    frequency_scale(model::Potential)

Frequency scale associated to the potential `model`.
"""
function frequency_scale(model::Potential)
    # ... [model specific implementation] ...
end

"""
    energy_scale(model::Potential)

Energy scale associated to the potential `model`.
"""
function energy_scale(model::Potential)
    # ... [model specific implementation] ...
end

"""
    momentum_scale(model::Potential)

Angular momentum scale associated to the potential `model`.
"""
function momentum_scale(model::Potential)
    # ... [model specific implementation] ...
end

"""
    radial_scale(model::Potential)

Radial scale associated to the potential `model`.
"""
function radial_scale(model::Potential)
    # ... [model specific implementation] ...
end

#####################################
# Include potentials
#####################################
include("isochrone.jl")
include("plummer.jl")
include("mestelzang.jl")