#####################################
# Abstract types
#####################################
abstract type OrbitalDistribution end

abstract type SphericalDistribution <: OrbitalDistribution end
abstract type EnergySphericalDistribution <: SphericalDistribution end
abstract type ActionsSphericalDistribution <: SphericalDistribution end

abstract type RazorThinDistribution <: OrbitalDistribution end
abstract type EnergyRazorDistribution <: RazorThinDistribution end
abstract type ActionsRazorDistribution <: RazorThinDistribution end


#####################################
# Include specific distribution functions
#####################################
include("isochrone.jl")
include("plummer.jl")
include("mestelzang.jl")
include("miyamoto.jl")