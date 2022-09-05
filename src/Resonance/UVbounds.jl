"""

"""

"""Getϖ(ω₀,n₁,n₂,dψ/dr,d²ψ/dr²[,Ω₀,rmin,rmax])
translate a complex frequency into a rescale frequency.
maps ``\\omega \\to [-1,1]``

Fouvry & Prunet B3

@ASSUMPTION:
    - ω is dimensionless, that is, rescaled by Ω₀ already.
"""
function Getϖ(ω::Complex{Float64},
                   n1::Int64,n2::Int64,
                   dψ::Function,d2ψ::Function;
                   Ω₀::Float64=1.,
                   rmin::Float64=1.e-10,
                   rmax::Float64=1000.)

    ωmin, ωmax = Findωminωmax(n1,n2,dψ,d2ψ;Ω₀=Ω₀,rmin=rmin,rmax=rmax)

    return Getϖ(ω,ωmin,ωmax)
end

"""
ϖ version with ωmin, ωmax

"""
function Getϖ(ω::Complex{Float64},
              ωmin::Float64,ωmax::Float64)

    return (2.0*ω - ωmax - ωmin)/(ωmax - ωmin)
end



########################################################################
#
# (u,v) mapping : ω boundaries (at given n1, n2)
#
########################################################################

"""Findωminωmax(n₁,n₂,dψ,d2ψ[,rmax,Ω₀])
for a given resonance, find the maximum frequencies

@ASSUMPTION:
    - Frenquency domain truncated at αmin and αmax
"""
function Findωminωmax(n1::Int64,n2::Int64,
                      dψ::Function,
                      d2ψ::Function;
                      Ω₀::Float64=1.,
                      rmin::Float64=1.0e-8,
                      rmax::Float64=1.0e5)

    # define the function to extremise
    ωncirc(x) = n1*Ω1circular(dψ,d2ψ,x)/Ω₀ + n2*Ω2circular(dψ,d2ψ,x)/Ω₀

    αmin, αmax = αminmax(dψ,d2ψ,rmin,rmax,Ω₀=Ω₀)
    xext = ExtremiseFunction(ωncirc,rmin,rmax)

    # The extreme values of n.Ω is either :
    #   - on the radial line, at α = αmin or αmax
    #   - along the circular velocity (extreme α included)
    ωmin,ωmax = extrema([ωncirc(xext), ωncirc(rmin), ωncirc(rmax), (n1+0.5*n2)*αmin, (n1+0.5*n2)*αmax])
    return ωmin,ωmax
end

"""αminmax(dψ,d2ψ,rmin,rmax[, Ω₀])
maximal and minimal considered radial frequencies (rescaled)

@ASSUMPTION:
    - Ω1circular is a decreasing function of radius
"""
function αminmax(dψ::Function,
                 d2ψ::Function,
                 rmin::Float64,
                 rmax::Float64;
                 Ω₀::Float64=1.)

    @assert rmin < rmax "rmin >= rmax in αminmax function"
    # Assumption :
    # Ω1circular is a decreasing function of radius
    return Ω1circular(dψ,d2ψ,rmax)/Ω₀, Ω1circular(dψ,d2ψ,rmin)/Ω₀
end

########################################################################
#
# (u,v) mapping : v boundary (at given u, n1, n2)
#
########################################################################



"""FindVminVmax(u,ωmin,ωmax,n₁,n₂,vbound,βc)
for a given resonance, at a specific value of u, find the v coordinate boundaries.

@IMPROVE, put in guards for the edges in βC
@ASSUMPTION:
    - rmin, rmax are the same used for ωmin, ωmax, αmin and αmax computation
"""
function FindVminVmax(u::Float64,
                      n1::Int64,n2::Int64,
                      dψ::Function,
                      d2ψ::Function,
                      ωmin::Float64,ωmax::Float64,
                      αmin::Float64,αmax::Float64,
                      βc::Function;
                      Ω₀::Float64=1.,
                      rmin::Float64=1.0e-8,
                      rmax::Float64=1.0e5)


    # ωn(u) : value of the resonance line
    hval = HUFunc(u,ωmin,ωmax)

    if (n2==0)
        #####
        # (B9) Fouvry & Prunet : 1st inequality
        #####
        vmin = 0.5

        #####
        # (B9) Fouvry & Prunet : 2nd inequality
        #####
        # put in guards for the edges
        αvmax = min(αmax,max(hval/n1,αmin))

        vmax = βc(αvmax)
    else
        #####
        # (B10) Fouvry & Prunet : 1st & 2nd inequalities
        #####
        vmin = αmin
        vmax = αmax

        #####
        # (B10) Fouvry & Prunet : 3rd inequality
        #####
        kappa = n1+0.5*n2
        if (n2*hval > 0.)
            if (kappa*hval > 0.)
                vmax = min(vmax, hval/kappa) # Updating vmax
            end
        elseif (n2*hval < 0.)
            if     (kappa*hval > 0.)
                vmin = max(vmin, hval/kappa) # Updating vmin
            end
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
            vbound = FindVbound(n1,n2,dψ,d2ψ,Ω₀=Ω₀,rmin=rmin,rmax=rmax)

            # Extreme boundary to look for vbound
            # @WARNING arbitrary fixed constant
            locrmin, locrmax = 1.e-6, 1.e5

            if (vbound != αmin) && (vbound != αmax)
                branch = 2
            elseif (rmin > locrmin) || (rmax > locrmax)
                # If vboung not in the asked boundary
                # verify that it should indeed not exist
                vbound = FindVbound(n1,n2,dψ,d2ψ,Ω₀=Ω₀,rmin=locrmin,rmax=locrmax)

                if (vbound != Ω1circular(dψ,d2ψ,locrmin)/Ω₀) && (vbound != Ω1circular(dψ,d2ψ,locrmin)/Ω₀)
                    branch = 2
                else
                    branch = 1
                end
            else
                branch = 1
            end
        end

        # Constraint equation
        rootequation(v::Float64)::Float64 = hval - n1*v - n2*v*βc(v)

        if branch == 1 # ωncirc(x) is monotonic
            if rootequation(αmin)*rootequation(αmax) < 0.
                vlim = bisection(rootequation,αmin,αmax)
                if hval*n2 > 0.
                    vmin = max(vmin,vlim)
                elseif hval*n2 < 0.
                    vmax = min(vmax,vlim)
                end
            end
        elseif branch == 2 # ωncirc(x) is not monotonic
            # Search crossing in [αmin,vbound]
            if rootequation(αmin)*rootequation(vbound) < 0.
                vmin2 = bisection(rootequation,αmin,vbound)
                vmin = max(vmin,vmin2)
            end
            # Search crossing in [vbound,αmax]
            if rootequation(vbound)*rootequation(αmax) < 0.
                vmax2 = bisection(rootequation,vbound,αmax)
                vmax = min(vmax,vmax2)
            end
        end
    end

    return vmin, vmax
end

"""HUFunc(u,ωmin,ωmax)
return h_n(u) = ω_n(u), a helper quantity
Fouvry & Prunet B8
"""
function HUFunc(u::Float64,ωmin::Float64,ωmax::Float64)
    return 0.5*(ωmax+ωmin + u*(ωmax-ωmin))
end


"""FindVbound(n₁,n₂,dψ,d2ψ[,rmax,Ω₀])
find any valie non- 0 or 1 v value at u=-1 or u=1
"""
function FindVbound(n1::Int64,n2::Int64,
                    dψ::Function,
                    d2ψ::Function;
                    Ω₀::Float64=1.,
                    rmin::Float64=1.0e-8,
                    rmax::Float64=1.0e5)

    # define the function to extremise
    ωncirc(x) = n1*Ω1circular(dψ,d2ψ,x) + n2*Ω2circular(dψ,x)

    xext = ExtremiseFunction(ωncirc,rmin,rmax)

    return Ω1circular(dψ,d2ψ,xext)/Ω₀
end
