# structure to hold all parameters
include("parameters.jl")

# utils: simple function extremisation, input/output, basic orbit coordinate transformations, simple integrators
include("IO.jl")
include("edges.jl")
include("extremise.jl")
include("integrate.jl")
include("differentiate.jl")
include("interpolate.jl")
include("backward.jl")
include("spheroidalcoordinates.jl")