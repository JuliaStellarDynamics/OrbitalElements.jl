"""
The plummer potential definitions

  plummer.jl is an example of how to implement a simple new function and derivatives.
  a frequency scaling is also required.
"""

#####################################
# Plummer structures
#####################################
abstract type PlummerPotential <: TwoIntegralCentralCorePotential end
struct NumericalPlummer <: PlummerPotential
    G::Float64      # Gravitational constant
    M::Float64      # Total mass
    bc::Float64     # Characteristic radius
end
struct SemiAnalyticPlummer <: PlummerPotential
    G::Float64      # Gravitational constant
    M::Float64      # Total mass
    bc::Float64     # Characteristic radius
end

"""
    NumericalPlummer([, G, M, bc])

This plummer model will use the default numerical computations
"""
function NumericalPlummer(;G::Float64=1.,M::Float64=1.,bc::Float64=1.)
    # Check for positive mass and radius
    if M<0; throw(DomainError(M, "Negative mass")); end
    if bc≤0; throw(DomainError(bc, "Negative characteristic radius")); end

    return NumericalPlummer(G,M,bc)
end

"""
    SemiAnalyticPlummer([, G, M, bc])

This Plummer model will use the semi-analytic mappings with an anomaly 
such that the frequency integrand [`Θ(u, a, e, ...)`](@ref) is known.
The frequencies are still computed through numerical integration of this integrand.
"""
function SemiAnalyticPlummer(;G::Float64=1.,M::Float64=1.,bc::Float64=1.)
    # Check for positive mass and radius
    if M<0; throw(DomainError(M, "Negative mass")); end
    if bc≤0; throw(DomainError(bc, "Negative characteristic radius")); end

    return SemiAnalyticPlummer(G,M,bc)
end

#####################################
# Potential methods for Plummer
#####################################
function ψ(r::Float64,model::PlummerPotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc # Rescaled radius
    scale = model.G*model.M/model.bc # Scale
    return - scale / sqrt(1.0+x^2)
end

function dψ(r::Float64,model::PlummerPotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end
    
    x = r/model.bc
    scale = model.G*model.M/(model.bc^2)
    # Stable version at infinity (not stable in 0.)
    if x > 1.e5
        return scale / (x^2 * (sqrt(1.0+x^(-2)))^3)
    end
    return scale * x / (sqrt(1.0+x^2))^3
end

function d2ψ(r::Float64,model::PlummerPotential)
    # Check for positive radius
    if r<0; throw(DomainError(r, "Negative radius")); end

    x = r/model.bc
    scale = model.G*model.M/(model.bc^3)
    # Stable version at infinity (not stable in 0.)
    if x > 1.e5
        return scale * ( 1.0 / (sqrt(1.0 + x^2))^3 
                        - 3.0 / (x^3 * (sqrt(1.0 + x^(-2)))^5))
    end
    return scale * (1.0 - 2.0*(x^2)) / (sqrt(1.0 + x^2))^5
end

#####################################
# Scales for Plummer
#####################################
function frequency_scale(model::PlummerPotential)
    return 2*sqrt(model.G*model.M/(model.bc^3))
end

"""
for Plummer, see Tep+ 22 (equation E2).
"""
function energy_scale(model::PlummerPotential)
    return -model.G*model.M/model.bc
end

"""
for Plummer, see Tep+ 22 (equation E2).
"""
function momentum_scale(model::PlummerPotential)
    return sqrt(model.G*model.M*model.bc)
end

function radial_scale(model::PlummerPotential)
    return model.bc
end