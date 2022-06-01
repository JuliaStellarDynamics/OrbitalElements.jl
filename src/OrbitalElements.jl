module OrbitalElements

# for zero-finding
using Optim
using Roots

# for interpolations
using Interpolations

# utils
include("Utils/Extremise.jl")
include("Utils/IO.jl")
include("Utils/OrbitDefinitions.jl")

# bring in the test potentials (not strictly needed)
include("Potentials/isochrone.jl")
include("Potentials/plummer.jl")

# enable energy and angular momentum computation (including expansions)
# first, set the radius where we switch to expansions
const ELTOLECC = 0.00001
include("Utils/ComputeEL.jl")

# bring in the circular orbit frequencies
include("CircularRadial/SpecialCases.jl")

# bring in the resonance mappings
include("Resonance/UVbounds.jl")
include("Resonance/ABtoUV.jl")

# the main wrapper for frequency and action calculations
include("Frequencies.jl")

end # module
