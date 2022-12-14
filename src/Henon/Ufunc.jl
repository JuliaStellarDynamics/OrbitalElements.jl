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

"""the henon anomaly increment
"""
function henonf(u::Float64)::Float64
    return u*(1.5 - 0.5*(u^2))
end

"""the derivative of the henon anomaly increment
"""
function henondf(u::Float64)::Float64
    return 1.5*(1.0 - u^(2))
end

"""the second derivative of the henon anomaly increment
"""
function henond2f(u::Float64)::Float64
    return -3u
end

"""the third derivative of the henon anomaly increment
"""
function henond3f(u::Float64)::Float64
    return -3.
end

"""the fourth derivative of the henon anomaly increment
"""
function henond4f(u::Float64)::Float64
    return 0.
end


########################################################################
#
# Radius as a function of u and derivatives in Hénon anomaly mapping
#
########################################################################

"""ru(u,a,e)
mapping from u->r in Henon variables

"""
function ru(u::Float64,a::Float64,e::Float64)::Float64
    return a*(1+e*henonf(u))
end


"""drdu(u,a,e)
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
the effective potential: note the relationship to Q
"""
function ψeff(ψ::Function,r::Float64,L::Float64)::Float64
    if L == 0.
        return ψ(r)
    else
        return ψ(r) + 0.5 * (L/r)^(2)
    end
end

"""
the derivative of the effective potential
"""
function dψeffdr(dψ::Function,r::Float64,L::Float64)::Float64
    if L == 0.
        return dψ(r)
    else
        return dψ(r) - (L)^(2) / (r^3)
    end
end

"""
the second derivative of the effective potential
"""
function d2ψeffdr2(d2ψ::Function,r::Float64,L::Float64)::Float64
    if L == 0.
        return d2ψ(r)
    else
        return d2ψ(r) + 3 * (L)^(2) / (r^4)
    end
end


########################################################################
#
# Θ(u) (frequencies integrand anomaly)
#
########################################################################

"""ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e[,EDGE,TOLECC,f,df,d2f,d3f,d4f])

Θ, the anomaly for computing orbit averages as a function of (a,e)

equivalent to Θ = (dr/du)(1/Vrad)

"""
function ΘAE(ψ::F0,
             dψ::F1,
             d2ψ::F2,
             d3ψ::F3,
             u::Float64,
             a::Float64,
             e::Float64,
             params::OrbitsParameters)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # use the expanded approximation
    if ((1-abs(u))<params.EDGE)
        return ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e,params)

    # if not close to the boundary, can calculate as normal
    else
        E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)

        r = ru(u,a,e)

        # this can somehow be negative: do we need an extra check?
        denomsq = 2*(E - ψeff(ψ,r,L))

        if (denomsq <= 0.0) || isnan(denomsq) || isinf(denomsq)
            # go back to the expansion -- or should we return 0.0?
            return ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e,params)
        else 
            return a * e * henondf(u) / sqrt(denomsq)
        end
    end
end




"""ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e[,EDGE,TOLECC,f,df,d2f,d3f,d4f])

ΘExpansion, the anomaly for computing orbit averages as a function of (a,e)

Used when u is sufficiently close to +1,-1

BIG ALLOCATIONS here from ELTOLECC not being specified.
Downside is that this guarantees allocations if TOLECC not specified.

"""
function ΘExpansionAE(ψ::F0,
                        dψ::F1,
                        d2ψ::F2,
                        d3ψ::F3,
                        u::Float64,
                        a::Float64,
                        e::Float64,
                        params::OrbitsParameters)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # which boundary are we close to?
    ul = (u > 0.) ? 1.0 : -1.0

    # compute the corresponding radius value
    rl = ru(ul,a,e)

    # compute energy and angular momentum from the potential (allow for expansions)
    L = LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)

    # compute the derivatives of the effective potential
    dψeffl, d2ψeffl = dψeffdr(dψ,rl,L), d2ψeffdr2(d2ψ,rl,L)

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
        return 0.0
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


"""ΘAEdade(ψ,dψ,d2ψ,d3ψ,u,a,e[,EDGE,TOLECC,f,df,d2f,d3f,d4f,da,de])

numerical differentiation of Θ w.r.t. semimajor axis and eccentricity
"""
function ΘAEdade(ψ::F0,
                     dψ::F1,
                     d2ψ::F2,
                     d3ψ::F3,
                     u::Float64,
                     a::Float64,
                     e::Float64,
                     params::OrbitsParameters)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    da, de = params.da, params.de

    # derivative w.r.t. semimajor axis: always safe
    thHa = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a+da,e,params)
    thLa = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a   ,e,params)
    dΘda = (thHa-thLa)/da

    # derivative w.r.t. semimajor axis: safe unless too close to e=1
    # in that case, reverse de
    if e>(1-de)
        de *= -1
    end

    thHe = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e+de,params)
    thLe = thLa
    dΘde = (thHe-thLe)/de

    return dΘda,dΘde
end