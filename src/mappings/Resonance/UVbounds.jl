"""

"""

"""
    Getϖ(ω,n1,n2,model,params)

translate a complex frequency into a rescaled frequency.
maps ``\\omega \\to [-1,1]``

Fouvry & Prunet B3

@ASSUMPTION:
    - ω is dimensionless, that is, rescaled by Ω₀ already.
"""
function Getϖ(ω::Number,
              n1::Int64,n2::Int64,
              model::CentralPotential,
              params::OrbitalParameters=OrbitalParameters())

    ωmin, ωmax = Findωminωmax(n1,n2,model,params)

    return Getϖ(ω,ωmin,ωmax)
end

"""
    Getϖ(ω,ωmin,ωmax)

ϖ version with ωmin, ωmax

"""
function Getϖ(ω::Number,
              ωmin::Float64,ωmax::Float64)

    return (2.0*ω - ωmax - ωmin)/(ωmax - ωmin)
end


########################################################################
#
# (u,v) mapping : ω boundaries (at given n1, n2)
#
########################################################################

"""
    Findωminωmax(n1,n2,model,params)
for a given resonance, find the maximum frequencies

@ASSUMPTION:
    - Frenquency domain truncated at αmin and αmax
"""
function Findωminωmax(n1::Int64,n2::Int64,
                      model::CentralPotential,
                      params::OrbitalParameters)::Tuple{Float64,Float64}

    Ω₀ = params.Ω₀
    rmin, rmax = params.rmin, params.rmax
    αmin, αmax = αminmax(model,rmin,rmax,Ω₀)
    
    # define the function to extremise
    ωncirc(x::Float64)::Float64 = n1*Ω1circular(model,x)/Ω₀ + n2*Ω2circular(model,x)/Ω₀

    # If rmax is infinite, bisection search on a bounded interval
    xext = ExtremiseFunctionNulCure(ωncirc,rmin,min(rmax,1.e8*params.rc))

    # The extreme values of n.Ω is either :
    #   - on the radial line, at α = αmin or αmax
    #   - along the circular velocity (extreme α included)
    ωmin = min(ωncirc(xext), ωncirc(rmin), ωncirc(rmax), (n1+0.5*n2)*αmin, (n1+0.5*n2)*αmax)
    ωmax = max(ωncirc(xext), ωncirc(rmin), ωncirc(rmax), (n1+0.5*n2)*αmin, (n1+0.5*n2)*αmax)

    return ωmin, ωmax
end

"""
    αminmax(model,rmin,rmax,Ω₀)

maximal and minimal considered radial frequencies (rescaled)

@ASSUMPTION:
    - Ω1circular is a decreasing function of radius
"""
function αminmax(model::CentralPotential,
                 rmin::Float64,rmax::Float64,
                 Ω₀::Float64)::Tuple{Float64,Float64}

    @assert rmin < rmax "rmin >= rmax in αminmax function"
    # Assumption :
    # Ω1circular is a decreasing function of radius
    return Ω1circular(model,rmax)/Ω₀, Ω1circular(model,rmin)/Ω₀
end



########################################################################
#
# (u,v) mapping : v boundary (at given u, n1, n2)
#
########################################################################

"""
    FindVminVmax(u,n1,n2,model,ωmin,ωmax,params)

for a given resonance, at a specific value of u, find the v coordinate boundaries.

@ASSUMPTION:
    - rmin, rmax are the same used for ωmin, ωmax, αmin and αmax computation
"""
function FindVminVmax(u::Float64,
                      n1::Int64,n2::Int64,
                      model::CentralPotential,
                      ωmin::Float64,ωmax::Float64,
                      params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64}


    Ω₀ = params.Ω₀
    rc = params.rc
    rmin, rmax = params.rmin, params.rmax
    αmin, αmax = αminmax(model,rmin,rmax,Ω₀)

    # ωn(u) : value of the resonance line
    hval = HUFunc(u,ωmin,ωmax)

    # βcircular as a function of αcircular
    βc(αc::Float64)::Float64 = βcirc(αc,model,params)

    if (n2==0) # v = β
        #####
        # (B9) Fouvry & Prunet : 1st inequality
        #####
        vmin = 0.5

        #####
        # (B9) Fouvry & Prunet : 2nd inequality
        #####
        vmax = βc(hval/n1)
    else # v = α

        #####
        # (B10) Fouvry & Prunet : 1st & 2nd inequalities
        #####
        vmin = αmin
        vmax = αmax

        #####
        # Constraint for hval = 0.
        # Line of constant β (or v = 0.)
        #####
        # β = - n1/n2 ≥ 1/2
        if (hval == 0.) && (-n1/n2 < 0.5)
            vmin = 0.
            vmax = 0.
        end
        # β = - n1/n2 ≤ βc(v)
        # Hypothesis : βc is a decreasing function of α (=v)
        if (hval == 0.) && (-n1/n2 >= βc(αmin))
            vmin = αmin
            vmax = αmin
        end

        #####
        # (B10) Fouvry & Prunet : 3rd inequality
        #####
        # Hiting radial boundary β = 1/2
        # hval = (n1+n2/2)v
        radon = n1+0.5*n2 # Radial orbit equivalent n
        if (n2*hval > 0.) && (radon*hval > 0.)
            vmax = min(vmax, hval/radon) # Updating vmax
        elseif (n2*hval < 0.) && (radon*hval > 0.)
            vmin = max(vmin, hval/radon) # Updating vmin
        end

        #####
        # (B10) Fouvry & Prunet : 4th inequality
        #####
        # branch == 1: ωncirc(x) is monotonic
        # branch == 2: ωncirc(x) is not monotonic
        if n1*n2 > 0 # Do not search for vbound
            branch = 1
        else
            # First look for vbound in the asked boundary
            vbound = FindVbound(n1,n2,model,Ω₀,rc,rmin,rmax)

            # Extreme boundary to look for vbound
            locrmin, locrmax = 0., Inf

            if (vbound != αmin) && (vbound != αmax)
                branch = 2
            elseif (rmin > locrmin) || (rmax < locrmax)
                # If vbound not in the asked boundary
                # verify that it should indeed not exist
                locrmin, locrmax = min(rmin,locrmin), max(rmax,locrmax)
                vbound = FindVbound(n1,n2,model,Ω₀,rc,locrmin,locrmax)

                branch = ((vbound != Ω1circular(model,locrmin)/Ω₀) && (vbound != Ω1circular(model,locrmax)/Ω₀)) ? 2 : 1
            else
                branch = 1
            end
        end

        # Constraint equation
        rootequation(v::Float64)::Float64 = hval - n1*v - n2*v*βc(v)

        if branch == 1 # ωncirc(x) is monotonic
            vlim = try bisection(rootequation,αmin,αmax) catch; -1. end
            if (vlim != -1.) && (hval*n2 > 0.)
                vmin = max(vmin,vlim)
            elseif (vlim != -1.) && (hval*n2 < 0.)
                vmax = min(vmax,vlim)
            end
        elseif branch == 2 # ωncirc(x) is not monotonic
            # Search crossing in [αmin,vbound]
            vmin2 = try bisection(rootequation,αmin,vbound) catch; -1. end
            if (vmin2 != -1.)
                vmin = max(vmin,vmin2)
            end
            # Search crossing in [vbound,αmax]
            vmax2 = try bisection(rootequation,vbound,αmax) catch; -1. end
            if (vmax2 != -1.)
                vmax = min(vmax,vmax2)
            end
        end
    end

    return vmin, vmax
end


"""
    HUFunc(u,ωmin,ωmax)

return h_n(u) = ω_n(u), a helper quantity
Fouvry & Prunet B8
"""
function HUFunc(u::Float64,ωmin::Float64,ωmax::Float64)::Float64

    if u == 1.
        return ωmax
    elseif u == -1.
        return ωmin
    end 

    return 0.5*(ωmax+ωmin + u*(ωmax-ωmin))
end


"""
    FindVbound(n1,n2,model,Ω₀,rmin,rmax)

find any valid non- 0 or 1 v value at u=-1 or u=1
"""
function FindVbound(n1::Int64,n2::Int64,
                    model::CentralPotential,
                    Ω₀::Float64,rc::Float64,
                    rmin::Float64,rmax::Float64)::Float64

    # define the function to extremise
    ωncirc(x) = n1*Ω1circular(model,x) + n2*Ω2circular(model,x)

    # If rmax is infinite, bisection search on a bounded interval
    locrmax = min(rmax,1.e8*rc)
    xext = ExtremiseFunctionNulCure(ωncirc,rmin,locrmax)

    # If the extremum is reached at the imposed maximal boundary
    # Use the true rmax (not the artificial 1.e8*rc, which is here to handle Inf)
    if (xext == locrmax)
        xext = rmax
    end

    return Ω1circular(model,xext)/Ω₀
end
