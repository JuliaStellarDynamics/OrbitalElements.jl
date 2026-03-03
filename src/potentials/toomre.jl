
#####################################
#
# Classic Toomre potential
#
#####################################

#####################################
# Toomre structure
#####################################
"""
Toomre potential structure
"""
struct ToomrePotential <:  TwoIntegralCentralCorePotential  end
    G::Float64      # Gravitational constant
    M::Float64      # Total mass
    bc::Float64     # Characteristic radius
end

"""
    ToomrePotential([, R0, V0])

Create a Toomre potential structure. 
"""
function ToomrePotential(;G::Float64=1.,M::Float64=1.,bc::Float64=1.)
    return ToomrePotential(G,M,bc)
end

#####################################
# Potential methods for Toomre
#####################################
function ψ(r::Float64,model::ToomrePotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end
    
    x = r/model.bc
    scale = model.G * model.M / model.bc
    return -scale / sqrt(1.0 + x^2)
end

function dψ(r::Float64,model::ToomrePotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    scale = model.G * model.M / (model.bc)^2
    # Stable version at infinity (not stable in 0.)
    if x > 1.e5
        return scale / (x^2 * (sqrt(1.0+x^(-2)))^3)
    end
    return scale * x / (sqrt(1.0+x^2))^3
end

function d2ψ(r::Float64,model::ToomrePotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    scale = model.G * model.M / (model.bc)^3
    # Stable version at infinity (not stable in 0.)
    if x > 1.e5
        return scale * ( 1.0 / (sqrt(1.0 + x^2))^3 
                        - 3.0 / (x^3 * (sqrt(1.0 + x^(-2)))^5))
    end
    return scale * (1.0 - 2.0*(x^2)) / (sqrt(1.0 + x^2))^5
end

#####################################
# Scales for Toomre
#####################################
function frequency_scale(model::ToomrePotential)
    return sqrt(model.G*model.M/(model.bc^3))
end

function radial_scale(model::ToomrePotential)
    return model.bc
end

function energy_scale(model::ToomrePotential)
    return -model.G*model.M/model.bc
end

function momentum_scale(model::ToomrePotential)
    return sqrt(model.G*model.M*model.bc)
end