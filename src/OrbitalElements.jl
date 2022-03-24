module OrbitalElements

# for zero-finding
using Optim

# for interpolations
using Interpolations


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

# anything to IO?
include("IO.jl")

end # module
