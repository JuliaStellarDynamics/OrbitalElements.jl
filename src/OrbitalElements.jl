module OrbitalElements

#####################################
# Dependencies
#####################################
using HDF5                      # I/O
using HypergeometricFunctions   # potentials/distributions
using SpecialFunctions          # potentials/distributions

#####################################
# Exports
#####################################
# potentials
export Potential, CentralPotential, CentralCorePotential, CentralCuspPotential
export AnalyticIsochrone, NumericalIsochrone
export PlummerPotential
export MestelPotential, TaperedMestel
export ψ, dψ, d2ψ, energy_scale, momentum_scale, frequency_scale, radial_scale
# mappings (naive)
export rpra_from_ae, ae_from_rpra
# mappings (forward)
export EL_from_ae, actions_from_ae
export αβ_from_ae, frequencies_from_ae
# mappings (jacobians)
export ae_to_rpra_jacobian, rpra_to_ae_jacobian
# mappings (backward)
export ae_from_EL, ae_from_actions
export ae_from_αβ, ae_from_frequencies
# mappings (resonance)
export uv_from_αβ, αβ_from_ae
export Resonance
# methods
export OrbitalParameters

#####################################
# Includes
#####################################
# potentials models
include("potentials/potentials.jl")
# Include distribution functions
include("distributions/distributions.jl")
# Include mappings
include("mappings/mappings.jl")

end # module
