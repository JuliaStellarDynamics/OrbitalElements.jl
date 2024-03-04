
#####################################
# Generic functions (not methods)
#####################################
# Get frequencies
function frequencies_from_αβ(α,β,model) end
function frequencies_from_actions(J,L,model,params) end
function frequencies_from_ae(a,e,model,params) end
function frequencies_from_EL(E,L,model,params) end
function frequencies_from_rpra(rp,ra,model,params) end

# Get frequency ratios (α,β)
function αβ_from_actions(J,L,model,params) end
function αβ_from_ae(a,e,model,params) end
function αβ_from_EL(E,L,model,params) end
function αβ_from_frequencies(Ω₁,Ω₂,model) end
function αβ_from_rpra(rp,ra,model,params) end

# Get actions
function actions_from_αβ(α,β,model,params) end
function actions_from_ae(a,e,model,params) end
function actions_from_EL(E,L,model,params) end
function actions_from_frequencies(Ω₁,Ω₂,model,params) end
function actions_from_rpra(rp,ra,model,params) end

# Get energy and angular momentum
function EL_from_αβ(α,β,model,params) end
function EL_from_actions(J,L,model,params) end
function EL_from_ae(a,e,model,params) end
function EL_from_frequencies(Ω₁,Ω₂,model,params) end
function EL_from_rpra(rp,ra,model,params) end

# Get semi-major axis and eccentricity
function ae_from_αβ(α,β,model,params) end
function ae_from_actions(J,L,model,params) end
function ae_from_EL(E,L,model,params) end
function ae_from_frequencies(Ω₁,Ω₂,model,params) end
function ae_from_rpra(rp,ra) end

# Get peri and apocentre
function rpra_from_αβ(α,β,model,params) end
function rpra_from_actions(J,L,model,params) end
function rpra_from_EL(E,L,model,params) end
function rpra_from_frequencies(Ω₁,Ω₂,model,params) end
function rpra_from_ae(a,e) end



#####################################
# Includes
#####################################
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