# enable energy and angular momentum computation (including expansions)
include("ComputeEL.jl")

# bring in the circular orbit frequencies
include("Circular/CircularFrequencies.jl")


# bring in the resonance mappings
include("Resonance/UVbounds.jl")
include("Resonance/ResonanceLines.jl")
include("Resonance/ABtoUV.jl")


# the main wrapper for frequency and action calculations
include("Frequencies.jl")

# Add the analytic mappings
# @IMPROVE right now, HAS to be after general mapping to override
# the default methods
include("analytic.jl")