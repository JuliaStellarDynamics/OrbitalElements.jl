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
#####
# Potentials
#####
# potentials structures
export Potential, CentralPotential, CentralCorePotential, CentralCuspPotential
export IsochronePotential, AnalyticIsochrone, NumericalIsochrone
export PlummerPotential, NumericalPlummer, SemiAnalyticPlummer
export MestelPotential, TaperedMestel
# potential functions
export ψ, dψ, d2ψ
export energy_scale, momentum_scale, frequency_scale, radial_scale

#####
# Mappings
#####
# mappings (rp, ra)
export ae_from_rpra, rpra_from_ae
export ae_to_rpra_jacobian, rpra_to_ae_jacobian
# mappings (E, L)
export ae_from_EL, EL_from_ae 
export ae_to_EL_jacobian
# mappings actions
export ae_from_actions, actions_from_ae
# @IMPROVE: ae_to_actions_jacobian missing !
# export ae_to_actions_jacobian
# mappings (α, β)
export αβ_from_ae, ae_from_αβ
export ae_to_αβ_jacobian
# mappings frequencies
export ae_from_frequencies, frequencies_from_ae
export ae_to_frequencies_jacobian
# mappings (α, β) ↔ frequencies
export αβ_from_frequencies, frequencies_from_αβ
export αβ_to_frequencies_jacobian, frequencies_to_αβ_jacobian

# Resonant mappings structures
export Resonance
# Resonant mappings functions
export αβ_from_uv, uv_from_αβ
export uv_to_αβ_jacobian

# methods
export OrbitalParameters

#####
# Distributions
#####

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
