
#=
The Hernquist potential definitions

  hernquist.jl is an example of how to implement a simple new function and derivatives.
  a frequency scaling that creates unity frequencies at the centre is also required.

=#

abstract type HernquistPotential <: TwoIntegralCentralCuspPotential end
struct NumericalHernquist <: HernquistPotential
    G::Float64      # Gravitational constant
    M::Float64      # Total mass
    bc::Float64     # Characteristic radius
end

"""
    NumericalHernquist([, G, M, bc])

This Hernquist model will use the default numerical computations
"""
function NumericalHernquist(;G::Float64=1.,M::Float64=1.,bc::Float64=1.)
    # Check for positive mass and radius
    if M<0; throw(DomainError(M, "Negative mass")); end
    if bc≤0; throw(DomainError(bc, "Negative characteristic radius")); end

    return NumericalHernquist(G,M,bc)
end

"""
the Hernquist potential
"""
function ψ(r::Float64,model::HernquistPotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    return -(model.G*model.M/model.bc) / (1.0+x)
end

"""
the Hernquist potential derivative
"""
function dψ(r::Float64,model::HernquistPotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    return (model.G*model.M/(model.bc^2)) / ((1.0+x)^2)
end


"""
the Hernquist potential second derivative
"""
function d2ψ(r::Float64,model::HernquistPotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    return -(2*model.G*model.M/(model.bc^3)) / ((1.0+x)^3)
end


"""
a frequency scale for the Hernquist potential

@WARNING this is not the central frequency
"""
function frequency_scale(model::HernquistPotential)

    return sqrt(model.G*model.M/(model.bc^3))
end


function energy_scale(model::HernquistPotential)
    return -model.G*model.M/model.bc
end

"""
"""
function momentum_scale(model::HernquistPotential)
    return sqrt(model.G*model.M*model.bc)
end

function radial_scale(model::HernquistPotential)
    return model.bc
end