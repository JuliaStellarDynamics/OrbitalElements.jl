#####################################
# Includes
#####################################
# utility functions
include("utils/utils.jl")
# mappings
include("peri_apocentre.jl")
include("energy_momentum.jl")
include("frequencies.jl")
# Actions have to come after frequencies for anomaly/radial velocity definitions
include("actions.jl")
include("resonant.jl")

# Add the analytic mappings
# @IMPROVE right now, HAS to be after general mapping to overwrite
# the default methods
include("analytic.jl")