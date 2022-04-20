"""

@ IMPROVE, tabulate wmin,wmax once for all resonances of interest
"""


"""find_wmin_wmax(n1,n2,dpotential,ddpotential[,rmax,Omega0])

for a given resonance, find the maximum frequencies
must have Omega1_circular, Omega2_circular defined (CircularRadial/CircularFrequencies.jl)

@IMPROVE: specific to the cored cluster with infinite extent

OrbitalElements.find_wmin_wmax(-3,4,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr)
"""
function find_wmin_wmax(n1::Int64,n2::Int64,
                        dpotential::Function,
                        ddpotential::Function,
                        rmax::Float64=1000.,
                        Omega0::Float64=1.)

    # define the function to extremise
    extreme(x) = n1*Omega1_circular(dpotential,ddpotential,x) + n2*Omega2_circular(dpotential,x)

    # hard-coded to 24 iterations on extremise_function, but this could be a parameter for hyper accuracy
    m = extremise_function(extreme,24,0.,rmax,false)

    # for the problem of a cored cluster with an infinite extent, the w minima and maxima are either
    # in the very centre
    # in the very outskirts
    # along the circular velocity
    w_min,w_max = extrema([extreme(1.e-10)/Omega0,0.0,extreme(m)/Omega0])
    return w_min,w_max
end

"""get_varpi(omega,n1,n2,dpotential,ddpotential[,rmax,Omega0])

translate a complex frequency into a rescale frequency.
maps ``\\omega \\to [-1,1]``

Fouvry & Prunet B3

"""
function get_varpi(omg::Complex{Float64},n1::Int64,n2::Int64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.,Omega0::Float64=1.)

    w_min,w_max = find_wmin_wmax(n1,n2,dpotential,ddpotential,rmax,Omega0)

    return (2omg - w_max - w_min)/(w_max - w_min)

end

"""hu(u,wmin,wmax)

return h, a helper quantity

Fouvry & Prunet B8

"""
function hu(u::Float64,wmin::Float64,wmax::Float64)
    return 0.5*(wmax+wmin + u*(wmax-wmin))
end

"""root_of_h_omega(u,wmin,wmax,n1,n2,vound,beta_c)

solve for roots of the h(u) equation
Fouvry & Prunet B10 (term 3)

"""
function root_of_h_omega(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64,vbound::Float64,beta_c)

    hval = hu(u,wmin,wmax)
    rootequation(x) = hval - n1*x - n2*x*beta_c[x]
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

"""constraint_three(u,wmin,wmax,n1,n2)

helper function for finding v bounds
"""
function constraint_three(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64)
    # greater than this
    hval = hu(u,wmin,wmax)
    return hval/(n2/2 + n1)
end

"""find_vmin_vmax(u,wmin,wmax,n1,n2,vbound,beta_c)

for a given resonance, at a specific value of u, find the v coordinate boundaries.

@IMPROVE, put in guards for the edges in beta_c
@IMPROVE, is there a way to specify the type on beta_c?

"""
function find_vmin_vmax(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64,vbound::Float64,beta_c)
    # this function works for n2 != 0

    if (n2==0)

        hval = hu(u,wmin,wmax)

        vmin = 0.5

        # put in guards for the very edges. SLOPPY
        vmax = beta_c(minimum([0.99999,maximum([hval/n1,0.00001])]))

    else

        r1,r2 = root_of_h_omega(u,wmin,wmax,n1,n2,vbound,beta_c)

        r3 = constraint_three(u,wmin,wmax,n1,n2)

        if (abs(n2)<=abs(n1)) | (abs(n2)>=abs(2n1))
            vmin = maximum([0.0,r1])
            vmax = minimum([r2,1.0,r3])
        else
            vmin = maximum([0.0,r1,r3])
            vmax = minimum([r2,1.0])
        end

        if isnan(vmax)
            vmax = 1.
        end

        if isnan(vmin)
            vmin = 0.
        end
    end

    return vmin,vmax
end
