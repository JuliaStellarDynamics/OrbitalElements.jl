#####################################
# Abstract types
#####################################
abstract type OrbitalDistribution end

abstract type SphericalDistribution <: OrbitalDistribution end
abstract type EnergySphericalDF <: SphericalDistribution end
abstract type ActionsSphericalDF <: SphericalDistribution end

abstract type RazorThinDistribution <: OrbitalDistribution end
abstract type EnergyRazorDF <: RazorThinDistribution end
abstract type ActionsRazorDF <: RazorThinDistribution end


#####################################
# Include specific distribution functions
#####################################
include("isochrone.jl")
include("plummer.jl")
include("mestelzang.jl")
include("miyamoto.jl")