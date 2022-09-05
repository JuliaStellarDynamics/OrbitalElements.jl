"""
wrap the special cases of interest

"""

# circular orbits use the epicyclic approximation
include("CircularFrequencies.jl")

# all radial frequencies are natually well-handled:
# but could become a special case if needed
#include("RadialFrequencies.jl")
