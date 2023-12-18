module OrbitalElements

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
# Abstract types
#####################################
# recommendations for various default parameters:
const DEFAULT_TOLECC = 0.01
const DEFAULT_TOLA   = 0.1
const DEFAULT_EDGE   = 0.01
const DEFAULT_NINT   = 32
const DEFAULT_RMIN   = 0.
const DEFAULT_RC     = 1.
const DEFAULT_RMAX   = Inf
const DEFAULT_TOL    = 1.e-12
const DEFAULT_DA     = 1.e-4
const DEFAULT_DE     = 1.e-4
const DEFAULT_ITERMAX= 100
const DEFAULT_Ω0     = 1.0

# structure to hold all parameters
include("Utils/ParameterStructure.jl")

# utils: simple function extremisation, input/output, basic orbit coordinate transformations, simple integrators
include("Utils/Extremise.jl")
include("Utils/IO.jl")
include("Utils/OrbitDefinitions.jl")
include("Utils/Integrators.jl")
include("Utils/DerivationInterpolation.jl")
include("Utils/EdgeHandle.jl")

# enable energy and angular momentum computation (including expansions)
include("Utils/ComputeEL.jl")

# bring in the circular orbit frequencies
include("Circular/CircularFrequencies.jl")

# bring in the resonance mappings
include("Resonance/UVbounds.jl")
include("Resonance/ResonanceLines.jl")
include("Resonance/ABtoUV.jl")

# the main wrapper for frequency and action calculations
include("Frequencies.jl")

# bring in the test potentials (not strictly needed)
include("Potentials/isochrone.jl")
include("Potentials/plummer.jl")
include("Potentials/mestelzang.jl")

# bring in the test distribution functions (not strictly needed)
include("DistributionFunctions/isochrone.jl")
include("DistributionFunctions/plummer.jl")
include("DistributionFunctions/mestelzang.jl")
include("DistributionFunctions/miyamoto.jl")
include("DistributionFunctions/isochrone_discs.jl")

end # module
