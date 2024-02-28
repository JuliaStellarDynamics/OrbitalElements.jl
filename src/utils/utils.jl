# structure to hold all parameters
include("ParameterStructure.jl")

# utils: simple function extremisation, input/output, basic orbit coordinate transformations, simple integrators
include("Extremise.jl")
include("IO.jl")
include("OrbitDefinitions.jl")
include("Integrators.jl")
include("DerivationInterpolation.jl")
include("EdgeHandle.jl")