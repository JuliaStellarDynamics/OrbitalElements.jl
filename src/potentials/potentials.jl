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
    ψ(model::CentralPotential, r)

Potential value for central potential `model` at radius `r`.
"""
function ψ(model::CentralPotential, r::Number)
    # ... [model specific implementation] ...
end

"""
    dψ(model::CentralPotential, r)

Potential derivative value for central potential `model` at radius `r`.
"""
function dψ(model::CentralPotential, r::Number)
    # ... [model specific implementation] ...
end

"""
    d2ψ(model::CentralPotential, r)

Potential second derivative for central potential `model` at radius `r`.
"""
function d2ψ(model::CentralPotential, r::Number)
    # ... [model specific implementation] ...
end

"""
    Ω₀(model::Potential)

Frequency scale associated to the potential `model`.
"""
function Ω₀(model::Potential)
    # ... [model specific implementation] ...
end

"""
    E₀(model::Potential)

Energy scale associated to the potential `model`.
"""
function E₀(model::Potential)
    # ... [model specific implementation] ...
end

"""
    L₀(model::Potential)

Angular momentum scale associated to the potential `model`.
"""
function L₀(model::Potential)
    # ... [model specific implementation] ...
end

#####################################
# Include potentials
#####################################
include("isochrone.jl")
include("plummer.jl")
include("mestelzang.jl")