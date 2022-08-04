"""

@IMPROVE, tabulate wmin,wmax once for all resonances of interest
@IMPROVE, all w->ω?
"""


"""find_wmin_wmax(n₁,n₂,dψ/dr,d²ψ/dr²[,rmax,Ω0])
for a given resonance, find the maximum frequencies
must have Omega1_circular, Omega2_circular defined (CircularRadial/CircularFrequencies.jl)

@IMPROVE: specific to the cored cluster with infinite extent
@IMPROVE: fix rmax limits to be expected

OrbitalElements.find_wmin_wmax(-3,4,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr)
"""
function find_wmin_wmax(n1::Int64,n2::Int64,
                        dpotential::Function,
                        ddpotential::Function,
                        rmax::Float64=1000.,
                        Ω₀::Float64=1.;
                        Ziter=24)

    # define the function to extremise
    extreme(x) = n1*Omega1_circular(dpotential,ddpotential,x) + n2*Omega2_circular(dpotential,x)

    m = ExtremiseFunction(extreme,Ziter,0.,rmax,verbose=false)

    # for the problem of a cored cluster with an infinite extent, the w minima and maxima are either
    # in the very centre
    # in the very outskirts
    # along the circular velocity
    w_min,w_max = extrema([extreme(1.e-10)/Ω₀,0.0,extreme(m)/Ω₀])

    return w_min,w_max
end

"""find_vbound(n₁,n₂,dψ/dr,d²ψ/dr²[,rmax,Ω₀])
find any valie non- 0 or 1 v value at u=-1 or u=1
"""
function find_vbound(n1::Int64,n2::Int64,
                     dpotential::Function,
                     ddpotential::Function,
                     rmax::Float64=1000.,
                     Ω₀::Float64=1.;
                     Ziter=24)

    # define the function to extremise
    extreme(x) = n1*Omega1_circular(dpotential,ddpotential,x) + n2*Omega2_circular(dpotential,x)

    m = ExtremiseFunction(extreme,Ziter,0.,rmax,verbose=false)

    vbound = Omega1_circular(dpotential,ddpotential,m)/Ω₀

    return vbound

end


"""get_varpi(omega₀,n₁,n₂,dψ/dr,d²ψ/dr²[,rmax,Ω₀])
translate a complex frequency into a rescale frequency.
maps ``\\omega \\to [-1,1]``

Fouvry & Prunet B3

omg needs to come in dimensionless, that is, rescaled by Ω₀ already.

This is exact.
"""
function get_varpi(omg::Complex{Float64},
                   n1::Int64,n2::Int64,
                   dpotential::Function,ddpotential::Function;
                   rmax::Float64=1000.,Ω₀::Float64=1.)

    w_min,w_max = find_wmin_wmax(n1,n2,dpotential,ddpotential,rmax,Ω₀)

    return (2.0*omg - w_max - w_min)/(w_max - w_min)

end

"""
varpi version with w_min, w_max

"""
function get_varpi(omg::Complex{Float64},
                   w_min::Float64,w_max::Float64)

    return (2.0*omg - w_max - w_min)/(w_max - w_min)

end


"""hu(u,wmin,wmax)
return h, a helper quantity
Fouvry & Prunet B8
"""
function hu(u::Float64,wmin::Float64,wmax::Float64)
    return 0.5*(wmax+wmin + u*(wmax-wmin))
end

"""RootOfHOmega(u,wmin,wmax,n₁,n₂,vbound,βc)
solve for roots of the h(u) equation
Fouvry & Prunet B10 (term 3)

βc must be a function

@IMPROVE: rootequation can be optimised for memory footprint, perhaps?
@IMPROVE: is there a more clever way to do the root-finding?
@IMPROVE: decide if 1.e-6 tolerance at the ends of root-finding bound are safe
"""
function RootOfHOmega(u::Float64,
                         wmin::Float64,wmax::Float64,
                         n1::Int64,n2::Int64,
                         vbound::Float64,
                         βc::Function)

    hval = hu(u,wmin,wmax)

    rootequation(x::Float64)::Float64 = hval - n1*x - n2*x*βc(x)

    # two roots to try: bounded by [0,v(u=1)] and [v(u=1),1]
    if rootequation(0+1.e-6)*rootequation(vbound) < 0
        r1 = fzero(rootequation, 0+1.e-6,vbound)
    else
        r1 = 0.0
    end

    if rootequation(1-1.e-6)*rootequation(vbound) < 0
        r2 = fzero(rootequation, vbound,1-1.e-6)
    else
        r2 = 1.0
    end
    #print(r1," ",r2)

    # edge curing for vbounds
    # from upper left to lower right cure
    if (vbound > 0.999999) & (u==-1.) & (hu(-1.,wmin,wmax)<0)
        r1 = 0.999999
        r2 = 1.0
    end

    # from lower left to upper right cure
    if (vbound > 0.999999) & (u==1.) & (hu(1.,wmin,wmax)>0)
        r1 = 0.999999
        r2 = 1.
    end

    # greater than r1, less than r2
    return r1,r2
end

"""ConstraintThree(u,wmin,wmax,n₁,n₂)
helper function for finding v bounds

Fouvry & Prunet (2022) B10, third constraint

@ATTENTION: can return -Inf in n1=1,n2=-2 (or n1=-1,n2=2) scenario
"""
function ConstraintThree(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64)
    # greater than this
    hval = hu(u,wmin,wmax)
    denom = (n2/2 + n1)

    # special cases:
    if denom==0
        if n2 < 0
            return 1.
        else
            return 1.
        end
    else
        return hval/denom
    end
end


"""FindVminVmax(u,wmin,wmax,n₁,n₂,vbound,βc)
for a given resonance, at a specific value of u, find the v coordinate boundaries.

@IMPROVE, put in guards for the edges in BetaC
@IMPROVE, decide if we want to have FindVminVmax wrap the vbound calculation?
"""
function find_vmin_vmax(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64,vbound::Float64,βc::Function)
    # this function works for n2 != 0

    if (n2==0)

        hval = hu(u,wmin,wmax)

        vmin = 0.5

        # put in guards for the very edges. SLOPPY
        vmax = βc(minimum([0.99999,maximum([hval/n1,0.00001])]))

    else

        r1,r2 = RootOfHOmega(u,wmin,wmax,n1,n2,vbound,βc)

        r3 = ConstraintThree(u,wmin,wmax,n1,n2)

        if (abs(n2)<=abs(n1)) | (abs(n2)>=abs(2n1))
            vmin = maximum([0.0,r1])
            vmax = minimum([r2,1.0,r3])
        else
            vmin = maximum([0.0,r1,r3])
            vmax = minimum([r2,1.0])
        end

        if (isnan(vmax) | isinf(vmax))
            vmax = 1.
        end

        if (isnan(vmin) | isinf(vmin))
            vmin = 0.
        end
    end

    # account for reversed limits in some particular cases -- perhaps figure out why sometime.
    return min(vmin,vmax),max(vmin,vmax)
end
