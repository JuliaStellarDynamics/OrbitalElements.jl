module OrbitalElements

# define some basic orbit transformations
include("OrbitDefinitions.jl")
export rpra_from_ae

# bring in the test potentials
include("Potentials/isochrone.jl")
export isochrone_Omega_1_2

include("Potentials/plummer.jl")
# decide on what to export...

# enable energy and angular momentum computation (including expansions)
include("ComputeEL.jl")
export E_from_rpra_pot,L_from_rpra_pot

# bring in the circular orbit frequencies
include("CircularRadial/CircularFrequencies.jl")
export Omega1_expansion,Omega2_expansion

# bring in the generic frequency solver
include("Henon/Frequencies.jl")
export compute_frequencies_henon

# anything to IO?
include("IO.jl")

end # module
