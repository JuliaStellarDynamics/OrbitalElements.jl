"""
Ufunc.jl

collection of functions that are used to compute quantities when integrating over an orbit angle.

Specific to Henon mapping.

"""

########################################################################
#
# Hénon anomaly function and derivatives
#
########################################################################

"""
    henonf(u)

the henon anomaly increment
"""
function henonf(u::Float64)::Float64
    return u*(1.5 - 0.5*(u^2))
end

"""
    henondf(u)

the derivative of the henon anomaly increment
"""
function henondf(u::Float64)::Float64
    return 1.5*(1.0 - u^(2))
end

"""
    henond2f(u)
    
the second derivative of the henon anomaly increment
"""
function henond2f(u::Float64)::Float64
    return -3u
end

"""
    henond3f(u)

the third derivative of the henon anomaly increment
"""
function henond3f(u::Float64)::Float64
    return -3.
end

"""
    henond4f(u)

the fourth derivative of the henon anomaly increment
"""
function henond4f(u::Float64)::Float64
    return 0.
end


########################################################################
#
# Radius as a function of u and derivatives in Hénon anomaly mapping
#
########################################################################

"""
    ru(u,a,e)

mapping from u->r in Henon variables
"""
function ru(u::Float64,a::Float64,e::Float64)::Float64
    return a*(1+e*henonf(u))
end


"""
    drdu(u,a,e)

mapping from u->r in Henon variables
"""
function drdu(u::Float64,a::Float64,e::Float64)::Float64
    return a*e*henondf(u)
end


########################################################################
#
# Effective potential and derivatives w.r.t. the radius
#
########################################################################

"""
    ψeff(ψ,r,L)

the effective potential: note the relationship to Q
"""
function ψeff(model::Potential,r::Float64,L::Float64)::Float64
    if L == 0.
        return ψ(model,r)
    end
        
    return ψ(model,r) + 0.5 * (L/r)^(2)
end

"""
    dψeffdr(dψ,r,L)

the derivative of the effective potential
"""
function dψeffdr(model::Potential,r::Float64,L::Float64)::Float64
    if L == 0.
        return dψ(model,r)
    end
        
    return dψ(model,r) - (L)^(2) / (r^3)
end

"""
    d2ψeffdr2(d2ψ,r,L)

the second derivative of the effective potential
"""
function d2ψeffdr2(model::Potential,r::Float64,L::Float64)::Float64
    if L == 0.
        return d2ψ(model,r)
    end
        
    return d2ψ(model,r) + 3 * (L)^(2) / (r^4)
end

########################################################################
#
# Radial velocity
#
########################################################################

"""
    Vrad(ψ,dψ,u,a,e,params)

radial velocity as a function of the orbital constants (a,e) and the anomaly u
"""
function Vrad(model::Potential,
              u::Float64,a::Float64,e::Float64,
              params::OrbitalParameters=OrbitalParameters())::Float64

    E, L = ELFromAE(model,a,e,params)

    r = ru(u,a,e)

    vrSQ = 2*(E - ψeff(model,r,L))

    if (vrSQ < 0.0) || isnan(vrSQ) || isinf(vrSQ)
        return 0.0
    end
        
    return sqrt(vrSQ)
end


########################################################################
#
# Θ(u) (frequencies integrand anomaly)
#
########################################################################

"""
    ΘAE(ψ,dψ,d2ψ,u,a,e,params)

Θ, the anomaly for computing orbit averages as a function of (a,e)
equivalent to Θ = (dr/du)(1/Vrad)
"""
function ΘAE(model::Potential,
             u::Float64,a::Float64,e::Float64,
             params::OrbitalParameters=OrbitalParameters())::Float64

    # use the expanded approximation
    # CAUTION: 1-(1-EDGE) < EDGE is true ...
    # To prevent this → EDGE - eps(Float64)
    if ((1.0 - abs(u)) < (params.EDGE - eps(Float64)))
        return ΘExpansionAE(model,u,a,e,params)
    end

    dr = drdu(u,a,e)

    # this can somehow be negative: do we need an extra check?
    vr = Vrad(model,u,a,e,params)

    if (vr == 0.0)
        # go back to the expansion -- or should we return 0.0?
        return 0.0
    end 
        
    return dr / vr
end




"""
    ΘExpansionAE(ψ,dψ,d2ψ,u,a,e,params)

ΘExpansion, the anomaly for computing orbit averages as a function of (a,e)

Used when u is sufficiently close to +1,-1
"""
function ΘExpansionAE(model::Potential,
                      u::Float64,a::Float64,e::Float64,
                      params::OrbitalParameters=OrbitalParameters())::Float64

    # which boundary are we close to?
    ul = (u > 0.) ? 1.0 : -1.0

    # compute the corresponding radius value
    rl = ru(ul,a,e)

    # compute energy and angular momentum from the potential (allow for expansions)
    L = LFromAE(model,a,e,params)

    # compute the derivatives of the effective potential
    dψeffl, d2ψeffl = dψeffdr(model,rl,L), d2ψeffdr2(model,rl,L)

    # compute the derivatives of the Henon f function
    d2fl, d3fl, d4fl = henond2f(ul), henond3f(ul), henond4f(ul)

    # define the prefactor
    pref = - ul * a * e

    # this denominator can be negative?
    combination = - a * e * dψeffl * d2fl

    # switch to safety: don't contribute anything at this point
    # In particular for radial orbits with ψ(r) = - Inf in r = 0
    # combination = - Inf
    if combination <= 0.
        # Switch to
        # extension of the function close to the border u ~ ul using the
        # linear interpolation between using the points ul ± tolu and ul ± 2*tolu
        u1 = ul * (1.0 - params.EDGE)
        u2 = ul * (1.0 - 2*params.EDGE)
        u3 = ul * (1.0 - 3*params.EDGE)

        Θ1 = ΘAE(model,u1,a,e,params)
        Θ2 = ΘAE(model,u2,a,e,params)
        Θ3 = ΘAE(model,u3,a,e,params)

        return Interpolation2ndOrder(u,u1,Θ1,u2,Θ2,u3,Θ3)
    end

    # if >0, sqrt is safe, proceed
    denom = sqrt(combination)

    zeroorder   = d2fl
    firstorder  = d3fl / 3.0
    secondorder = (3.0 * (dψeffl * d2fl * d4fl - a * e * d2ψeffl * (d2fl)^(3)) -  dψeffl * (d3fl)^(2)) / (24.0 * dψeffl * d2fl)

    return pref / denom * ( zeroorder + firstorder * (u - ul) + secondorder * (u - ul)^(2) )
end


########################################################################
#
# Θ(u) derivatives w.r.t. a and e
#
########################################################################

"""
    dΘAE(ψ,dψ,d2ψ,u,a,e,params)

numerical differentiation of Θ w.r.t. semimajor axis and eccentricity
"""
function dΘAE(model::Potential,
              u::Float64,a::Float64,e::Float64,
              params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64}

    # Function to differentiate
    fun(atemp::Float64,etemp::Float64) = ΘAE(model,u,atemp,etemp,params)
    # Perform differentiation
    _, ∂Θ∂a, ∂Θ∂e = NumericalDerivativeAE(fun,a,e,params)

    return ∂Θ∂a, ∂Θ∂e
end