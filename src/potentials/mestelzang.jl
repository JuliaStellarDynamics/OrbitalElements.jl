
#####################################
#
# Classic Mestel potential (cuspy)
#
#####################################

#####################################
# Mestel structure
#####################################
"""
Mestel potential structure
"""
struct MestelPotential <: TwoIntegralCentralCuspPotential 
    R0::Float64
    V0::Float64
end

"""
    MestelPotential([, R0, V0])

Create a Mestel potential structure with 'characteristic radius' `R0` and circular velocity `V0`. 

The characteristic radius just set the potential offset, fixing ``ψ(R0)=0``. 
"""
function MestelPotential(;R0::Float64=1.,V0::Float64=1.)
    return MestelPotential(R0,V0)
end

#####################################
# Potential methods for Mestel
#####################################
function ψ(r::Float64,model::MestelPotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end
    
    x = r/model.R0
    scale = model.V0^2
    return scale * log(x)
end

function dψ(r::Float64,model::MestelPotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.R0
    scale = model.V0^2 / model.R0
    return scale / x
end

function d2ψ(r::Float64,model::MestelPotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.R0
    scale = model.V0^2 / (model.R0^2)
    return  - scale / (x^2)
end

#####################################
# Scales for Mestel
#####################################
function frequency_scale(model::MestelPotential)
    return model.V0 / model.R0
end

function radial_scale(model::MestelPotential)
    return model.R0
end

# @IMPROVE define energy and momentum scales


#####################################
#
# Tapered Mestel potential (core)
#
#####################################

#####################################
# Tapered Mestel structure
#####################################
"""
Tapered Mestel potential structure
"""
struct TaperedMestel <: TwoIntegralCentralCorePotential 
    R0::Float64
    V0::Float64
    ε0::Float64
end

"""
    TaperedMestel([, R0, V0, ε0])

Create a tapered Mestel potential structure with 'characteristic radius' `R0` and circular velocity `V0` and taper length scale `ε0`. 
"""
function TaperedMestel(;R0::Float64=1.,V0::Float64=1.,ε0::Float64=1.e-5)
    return TaperedMestel(R0,V0,ε0)
end

#####################################
# Potential methods for tapered Mestel
#####################################
function ψ(r::Float64,model::TaperedMestel)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.R0
    scale = model.V0^2
    return scale * log(x^2 + model.ε0^2) / 2
end

function dψ(r::Float64,model::TaperedMestel)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.R0
    scale = model.V0^2 / model.R0
    # Stable version at infinity (not stable in 0.)
    if x > 1.e5
        return scale / (x * (1 + (model.ε0/x)^2))
    end
    return scale * x / (model.ε0^2 + x^2)
end

function d2ψ(r::Float64,model::TaperedMestel)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.R0
    scale = model.V0^2 / (model.R0^2)
    # Stable version at infinity (not stable in 0.)
    if x > 1.e5
        d2 = (model.ε0/x)^2
        return scale * (d2 - 1) / ((x * (1 + d2))^2)
    end
    return  scale * (model.ε0^2 - x^2) / ((model.ε0^2 + x^2)^2)
end

#####################################
# Scales for tapered Mestel
#####################################
function frequency_scale(model::TaperedMestel)
    return model.V0 / model.R0
end

function radial_scale(model::TaperedMestel)
    return model.R0
end

# @IMPROVE define energy and momentum scales