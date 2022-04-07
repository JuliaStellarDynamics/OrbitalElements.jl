module OrbitalElements

# for zero-finding
using Optim

# for interpolations
using Interpolations

# utils
include("Extremise.jl")

# anything to IO?
include("IO.jl")

# define some basic orbit transformations
include("OrbitDefinitions.jl")

# bring in the test potentials (not strictly needed)
include("Potentials/isochrone.jl")
include("Potentials/plummer.jl")

# enable energy and angular momentum computation (including expansions)
include("ComputeEL.jl")

# bring in the circular orbit frequencies
include("CircularRadial/SpecialCases.jl")

# bring in the generic frequency solver
include("Frequencies.jl")

# bring in the resonance mappings
include("Resonance/UVbounds.jl")


end # module
