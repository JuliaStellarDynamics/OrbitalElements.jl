module OrbitalElements

# for zero-finding
using Optim
using Roots

# for interpolations
using Interpolations

# utils: simple function extremisation, input/output, basic orbit coordinate transformations, simple integrators
include("Utils/Extremise.jl")
include("Utils/IO.jl")
include("Utils/OrbitDefinitions.jl")
include("Utils/Integrators.jl")

# bring in the test potentials (not strictly needed)
include("Potentials/isochrone.jl")
include("Potentials/plummer.jl")
include("Potentials/mestelzang.jl")
include("Potentials/kuzmintoomre.jl")

# bring in the test distribution functions (not strictly needed)
include("DistributionFunctions/isochrone.jl")
include("DistributionFunctions/plummer.jl")
include("DistributionFunctions/mestelzang.jl")
include("DistributionFunctions/miyamoto.jl")
include("DistributionFunctions/isochrone_discs.jl")

# enable energy and angular momentum computation (including expansions)
# first, set the radius where we switch to expansions
const ELTOLECC = 0.001
include("Utils/ComputeEL.jl")

# bring in the circular orbit frequencies
include("Circular/CircuarFrequencies.jl")

# bring in the resonance mappings
include("Resonance/UVbounds.jl")
include("Resonance/ABtoUV.jl")

# the main wrapper for frequency and action calculations
include("Frequencies.jl")

end # module
