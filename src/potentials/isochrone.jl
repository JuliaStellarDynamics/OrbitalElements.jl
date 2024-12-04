"""
The isochrone potential definitions

Isochrone is unique in that we can define almost all mappings analytically.
The AnalyticIsochrone potential model takes advantage of these. The analytic mappings are 
detailed in mappings/analytic.
Most equations can be found in Fouvry, Hamilton, Rozier & Pichon (2021), Appendix G
"""

#####################################
# Isochrone structures
#####################################

abstract type IsochronePotential <: TwoIntegralCentralCorePotential end
struct NumericalIsochrone <: IsochronePotential
    G::Float64      # Gravitational constant
    M::Float64      # Total mass
    bc::Float64     # Characteristic radius
end
struct AnalyticIsochrone <: IsochronePotential
    G::Float64      # Gravitational constant
    M::Float64      # Total mass
    bc::Float64     # Characteristic radius
end

"""
    NumericalIsochrone([, G, M, bc])

This isochrone model will use the default numerical computations
"""
function NumericalIsochrone(;G::Float64=1.,M::Float64=1.,bc::Float64=1.)
    # Check for positive mass and radius
    if M<0; throw(DomainError(M, "Negative mass")); end
    if bc≤0; throw(DomainError(bc, "Negative characteristic radius")); end

    return NumericalIsochrone(G,M,bc)
end

"""
    AnalyticIsochrone([, G, M, bc])

This isochrone model will use the full analytic mappings.
"""
function AnalyticIsochrone(;G::Float64=1.,M::Float64=1.,bc::Float64=1.)
    # Check for positive mass and radius
    if M<0; throw(DomainError(M, "Negative mass")); end
    if bc≤0; throw(DomainError(bc, "Negative characteristic radius")); end

    return AnalyticIsochrone(G,M,bc)
end

#####################################
# Potential methods for the Isochrone
#####################################
function ψ(r::Float64,model::IsochronePotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    scale = model.G*model.M/model.bc
    return - scale / (1.0+sqrt(1.0+x^2))
end

function dψ(r::Float64,model::IsochronePotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    scale = model.G*model.M/(model.bc^2)
    return scale / (sqrt(1.0+x^(-2))*(1.0+sqrt(1.0+x^2))^2)
end

function d2ψ(r::Float64,model::IsochronePotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    scale = model.G*model.M/(model.bc^3)
    sqxp = 1.0 + x^2
    invx = x^(-1)
    invsqxp = 1.0 + invx^2
    return scale * (- ((sqxp^(3/2))*(invx+sqrt(invsqxp))^2)^(-1)
                    - 2*(invsqxp*(1.0+sqrt(sqxp))^3)^(-1)
                    + (sqrt(sqxp)*(1.0+sqrt(sqxp))^2)^(-1))
end

#####################################
# Scales for the Isochrone
#####################################
function frequency_scale(model::IsochronePotential)
    return sqrt(model.G*model.M/(model.bc^3))
end

"""
for isochrone, see Fouvry+21 (appendix G).
"""
function energy_scale(model::IsochronePotential)
    return -model.G*model.M/model.bc
end

"""
for isochrone, see Fouvry+21 (appendix G)
"""
function momentum_scale(model::IsochronePotential)
    return sqrt(model.G*model.M*model.bc)
end

function radial_scale(model::IsochronePotential)
    return model.bc
end