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
export HernquistPotential, NumericalHernquist
export MestelPotential, TaperedMestel
export KuzminKutuzovPotential, AnalyticKuzminKutuzov
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
export EL_from_ae_derivatives
export ae_to_EL_jacobian
# mappings actions
export ae_from_actions, actions_from_ae
# @IMPROVE: ae_to_actions_jacobian missing !
# export ae_to_actions_jacobian
# mappings (α, β)
export αβ_from_ae, ae_from_αβ
export αβ_from_ae_derivatives
export ae_to_αβ_jacobian
# mappings frequencies
export ae_from_frequencies, frequencies_from_ae
export frequencies_from_ae_derivatives
export ae_to_frequencies_jacobian
# mappings (α, β) ↔ frequencies
export αβ_from_frequencies, frequencies_from_αβ
export αβ_from_frequencies_derivatives, frequencies_from_αβ_derivatives
export αβ_to_frequencies_jacobian, frequencies_to_αβ_jacobian
# Angle to radius and anomaly
export radius_from_anomaly
export angles_gradient

# Resonant mappings structures
export Resonance
# Resonant mappings functions
export resonance_number, frequency_extrema
export v_boundaries
export αβ_from_uv, uv_from_αβ
export uv_to_αβ_jacobian
export v_from_vp, vp_from_v
export v_from_vp_derivative
# Resonance lines helpers
export actions_resonance_line, actions_resonance_line!

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
# Include mappings
include("mappings/mappings.jl")

end # module
