module OrbitalElements

#####################################
# Dependencies
#####################################
using HDF5                      # I/O
using SpecialFunctions          # potentials/distributions
using HypergeometricFunctions   # potentials/distributions

#####################################
# Export
#####################################
# Potentials
export IsochronePotential, PlummerPotential, MestelPotential, TaperedMestel
export ψ, dψ, d2ψ, E₀, L₀, Ω₀
# Mappings (naive)
export RpRaFromAE, AEFromRpRa
# Mappings (forward)
export EFromAE, LFromAE, ELFromAE
export ComputeActionsAE, ComputeFrequenciesAE
# Mappings (backward)
export AEFromEL
export ComputeAEFromActions, ComputeAEFromFrequencies
# Mappings (resonance)
# Methods
export OrbitalParameters
#####################################
# Includes
#####################################
# Include utility functions
include("utils/utils.jl")
# Include potentials
include("potentials/potentials.jl")
# Include distribution functions
include("distributions/distributions.jl")
# Include mappings
include("mappings/mappings.jl")

end # module
