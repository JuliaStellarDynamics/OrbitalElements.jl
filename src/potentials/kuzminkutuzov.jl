"""
The Kuzmin-Kutuzov potential definitions

"""

#####################################
# Kuzmin-Kutuzov structures
#####################################

abstract type KuzminKutuzovPotential <: ThreeIntegralCentralCorePotential end

struct AnalyticKuzminKutuzov <: KuzminKutuzovPotential
    G::Float64      # Gravitational constant
    M::Float64      # Total mass
    a::Float64      # Characteristic radius
    c::Float64      # Characteristic radius

    ta::Float64     # Calculated from a and c
    tc::Float64     # Calculated from a and c

    alpha::Float64   # Calculated from a
    gamma::Float64   # Calculated from c
    Delta2::Float64 # Calculated from a and c
    Delta::Float64  # Calculated from Delta2
    sqrtDelta::Float64 # Calculated from Delta # these are part of a spheroidal coordinate system: a new abstract type?
end

"""
    AnalyticKuzminkutuzov([, G, M, a, c])

Creates an `AnalyticKuzminkutuzov` potential with the specified gravitational constant G,
mass M, and characteristic radii a and c. Checks for positive values of mass and radii.

The default value is the isotropic isochrone distribution.
"""
function AnalyticKuzminKutuzov(; G::Float64=1., M::Float64=1., a::Float64=0.5, c::Float64=0.5)
    # Check for positive mass and radius
    if M < 0; throw(DomainError(M, "Negative mass")); end
    if a <= 0; throw(DomainError(a, "Non-positive characteristic radius a")); end
    if c <= 0; throw(DomainError(c, "Non-positive characteristic radius c")); end

    # Calculate ta and tc based on a and c
    ta = a / (a + c)
    tc = c / (a + c)

    alpha = -a^2
    gamma = -c^2

    Delta2 = abs(gamma-alpha)
    Delta = sqrt(Delta2)
    sqrtDelta = sqrt(Delta)

    return AnalyticKuzminKutuzov(G, M, a, c, ta, tc, alpha, gamma, Delta2, Delta, sqrtDelta)
end



#####################################
# Potential methods for Kuzminkutuzov
#####################################

"""
    ψ(R::Float64, z::Float64, model::KuzminKutuzovPotential) -> Float64

Compute the gravitational potential ψ at a given cylindrical coordinate (R, z) for a Kuzmin-Kutuzov potential model.

# Arguments
- `R::Float64`: The radial distance from the z-axis. Must be non-negative.
- `z::Float64`: The height above the equatorial plane.
- `model::KuzminKutuzovPotential`: The Kuzmin-Kutuzov potential model containing parameters `a`, `c`, `G`, and `M`.

# Returns
- `Float64`: The gravitational potential at the specified coordinates.

# Throws
- `DomainError`: If `R` is negative.

# Notes
- This version not stable for R=0.,z=0.
"""
function ψ(R::Float64,z::Float64,model::KuzminKutuzovPotential)

    # Check for positive radius
    if R<0; throw(DomainError(R, "Negative radius")); end

    lambda,nu = lambda_nu_from_R_z(R,z,model.a,model.c)

    scale = model.G*model.M # Scale
    return - scale / (sqrt(lambda) + sqrt(nu)) 


end

function dψ(R::Float64,z::Float64,model::KuzminKutuzovPotential)

    # Check for positive radius
    if R<0; throw(DomainError(R, "Negative radius")); end

    lambda,nu = lambda_nu_from_R_z(R,z,model.a,model.c)


    scale = model.G*model.M

    return 1.
end

function d2ψ(R::Float64,z::Float64,model::KuzminKutuzovPotential)

    lambda,nu = lambda_nu_from_R_z(R,z,model.a,model.c)
    # Check for positive radius
    # Check for positive radius
    if lambda<0; throw(DomainError(lambda, "Negative radius")); end
    if nu<0; throw(DomainError(nu, "Negative radius")); end

    scale = model.G*model.M

    return 1.
end

#####################################
# Scales for Plummer
#####################################
function frequency_scale(model::KuzminKutuzovPotential)
    #@IMPROVE: check this is correct?
    return sqrt(model.G*model.M/((model.a+model.c)^3))
end

"""
for Kuzminkutuzov, see Tep+ 25 (equation 20).
"""
function energy_scale(model::KuzminKutuzovPotential)
    return -(model.a+model.c)/(model.G*model.M)
end

"""
for Kuzminkutuzov, see Tep+ 25 (equation 20).
"""
function momentum_scale(model::KuzminKutuzovPotential)
    return sqrt(model.G*model.M*(model.a+model.c))

end

function radial_scale(model::KuzminKutuzovPotential)
    return model.a
end

