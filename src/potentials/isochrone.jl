"""
The isochrone potential definitions

  isochrone.jl is unique in that we can define many of the quantities analytically, which makes for a useful testbed for empirical potentials.

  nearly all quantities for the isochrone potential can be computed analytically

"""

#####################################
# Isochrone structures
#####################################

"""
Isochrone potential structure
"""
struct IsochronePotential <: CentralCorePotential
    G::Float64      # Gravitational constant
    M::Float64      # Total mass
    bc::Float64     # Characteristic radius
end

"""
    IsochronePotential([, G, M, bc])

Create an Isochrone potential structure with characteristic radius `bc`, total mass `M` and gravitational constant `G`.
"""
function IsochronePotential(;G::Float64=1.,M::Float64=1.,bc::Float64=1.)
    # Check for positive mass and radius
    if M<0; throw(DomainError(M, "Negative mass")); end
    if bc≤0; throw(DomainError(bc, "Negative characteristic radius")); end

    return IsochronePotential(G,M,bc)
end

#####################################
# Potential methods for the Isochrone
#####################################
function ψ(model::IsochronePotential,r::Float64)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    scale = model.G*model.M/model.bc
    return - scale / (1.0+sqrt(1.0+x^2))
end

function dψ(model::IsochronePotential,r::Float64)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    scale = model.G*model.M/(model.bc^2)
    return scale / (sqrt(1.0+x^(-2))*(1.0+sqrt(1.0+x^2))^2)
end

function d2ψ(model::IsochronePotential,r::Float64)
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
function Ω₀(model::IsochronePotential)
    return sqrt(model.G*model.M/(model.bc^3))
end

"""
    E₀(model::IsochronePotential)

for isochrone, see Fouvry 21 (appendix G).
"""
function E₀(model::IsochronePotential)
    return -model.G*model.M/model.bc
end

"""
    L₀(model::IsochronePotential)

for isochrone, see Fouvry 21 (appendix G)
"""
function L₀(model::IsochronePotential)
    return sqrt(model.G*model.M*model.bc)
end