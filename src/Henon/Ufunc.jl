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
function henon_f(u::Float64)
    u*(3/2 - (u^2)/2)
end

"""the derivative of the henon anomaly increment
"""
function henon_df(u::Float64)
    1.5*(1.0 - u^(2))
end

"""the second derivative of the henon anomaly increment
"""
function henon_d2f(u::Float64)
    return -3u
end

"""the third derivative of the henon anomaly increment
"""
function henon_d3f(u::Float64)
    return -3.
end

"""the fourth derivative of the henon anomaly increment
"""
function henon_d4f(u::Float64)
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
function ru(u::Float64,a::Float64,e::Float64;fun::Function=henon_f)
    return a*(1+e*fun(u))
end

"""drdu(u,a,e)
mapping from u->r in Henon variables

"""
function drdu(u::Float64,a::Float64,e::Float64;dfun::Function=henon_df)
    return a*e*dfun(u)
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
        return ψ(r) + (1/2) * (L/r)^(2)
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
function ΘAE(ψ::Function,
             dψ::Function,
             d2ψ::Function,
             d3ψ::Function,
             u::Float64,
             a::Float64,
             e::Float64;
             EDGE::Float64=0.01,
             TOLECC::Float64=ELTOLECC,
             f::Function=henon_f,
             df::Function=henon_df,
             d2f::Function=henon_d2f,
             d3f::Function=henon_d3f,
             d4f::Function=henon_d4f)::Float64

    # use the expanded approximation
    if ((1-abs(u))<EDGE)
        return ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e;TOLECC=TOLECC,f=f,d2f=d2f,d3f=d3f,d4f=d4f)

    # if not close to the boundary, can calculate as normal
    else
        E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)

        r = ru(u,a,e;fun=f)

        # this can somehow be negative: do we need an extra check?
        denomsq = 2*(E - ψeff(ψ,r,L))

        if denomsq < 0.0
            # go back to the expansion -- or should we return 0.0?
            return ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e;TOLECC=TOLECC,f=f,d2f=d2f,d3f=d3f,d4f=d4f)
        end

        # do the standard return
        return a * e * df(u) / sqrt(denomsq)
    end
end

"""ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e[,EDGE,TOLECC,f,df,d2f,d3f,d4f])

ΘExpansion, the anomaly for computing orbit averages as a function of (a,e)

Used when u is sufficiently close to +1,-1

"""
function ΘExpansionAE(ψ::Function,
                          dψ::Function,
                          d2ψ::Function,
                          d3ψ::Function,
                          u::Float64,
                          a::Float64,
                          e::Float64;
                          TOLECC::Float64=ELTOLECC,
                          f::Function=henon_f,
                          df::Function=henon_df,
                          d2f::Function=henon_d2f,
                          d3f::Function=henon_d3f,
                          d4f::Function=henon_d4f)::Float64

    # which boundary are we close to?
    ul = (u > 0) ? 1.0 : -1.0

    # compute the corresponding radius value
    rl = ru(ul,a,e,fun=f)

    # compute energy and angular momentum from the potential (allow for expansions)
    E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)

    # compute the derivatives of the effective potential
    dψeffl, d2ψeffl = dψeffdr(dψ,rl,L), d2ψeffdr2(d2ψ,rl,L)

    # compute the derivatives of the Henon f function
    d2fl, d3fl, d4fl = d2f(ul), d3f(ul), d4f(ul)

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
    secondorder = (3.0 * (dψeffl * d2fl * d4fl - a * e * d2ψeffl * (d2fl)^(3)) -  dψeffl * (d3fl)^(2)) / (24 * dψeffl * d2fl)

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
function ΘAEdade(ψ::Function,
                     dψ::Function,
                     d2ψ::Function,
                     d3ψ::Function,
                     u::Float64,
                     a::Float64,
                     e::Float64;
                     EDGE::Float64=0.01,
                     TOLECC::Float64=ELTOLECC,
                     f::Function=henon_f,
                     df::Function=henon_df,
                     d2f::Function=henon_d2f,
                     d3f::Function=henon_d3f,
                     d4f::Function=henon_d4f,
                     da::Float64=1.0e-8,
                     de::Float64=1.0e-8)::Floa64

    # derivative w.r.t. semimajor axis: always safe
    thHa = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a+da,e,EDGE=EDGE,TOLECC=TOLECC,f=f,df=df,d2f=d2f,d3f=d3f,d4f=d4f)
    thLa = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a   ,e,EDGE=EDGE,TOLECC=TOLECC,f=f,df=df,d2f=d2f,d3f=d3f,d4f=d4f)
    dΘda = (thHa-thLa)/da

    # derivative w.r.t. semimajor axis: safe unless too close to e=1
    # in that case, reverse de
    if e>(1-de)
        de *= -1
    end

    thHe = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e+de,EDGE=EDGE,TOLECC=TOLECC,f=f,df=df,d2f=d2f,d3f=d3f,d4f=d4f)
    thLe = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e   ,EDGE=EDGE,TOLECC=TOLECC,f=f,df=df,d2f=d2f,d3f=d3f,d4f=d4f)
    dΘde = (thHe-thLe)/de

    return dΘda,dΘde
end
