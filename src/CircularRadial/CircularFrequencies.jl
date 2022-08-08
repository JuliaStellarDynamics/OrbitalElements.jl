#=
CircularFrequencies.jl

Special treatment for circular orbits
  for circular orbits, we can use exact relations with potential derivatives

@IMPROVE, add circular energy mapping

=#

"""Omega1_circular(dψ/dr,d²ψ/dr²,r)
radial frequency for circular orbits, from the epicyclic approximation

"""
function Omega1_circular(dpotential::Function,ddpotential::Function,r::Float64)
    return sqrt(ddpotential(r) + 3*dpotential(r)/r)
end

"""Omega2_circular(dψ/dr,r)
azimuthal frequency for circular orbits, from the epicyclic approximation
"""
function Omega2_circular(dpotential::Function,r::Float64)

    return sqrt(dpotential(r)/r)
end

"""Omega1_expansion(dψ/dr,d²ψ/dr²,r,rcirc)
first-order Taylor expansion for the radial frequency: at da
"""
function Omega1_expansion(dpotential::Function,ddpotential::Function,
                          r::Float64,rcirc::Float64,
                          dddpotential::Function=f(x)=x+1)

    omega1_circular = Omega1_circular(dpotential,ddpotential,r)
    dudr            = dpotential(rcirc)
    dduddr          = ddpotential(rcirc)
    h               = r-rcirc

    # did the user provide a third derivative?
    if dddpotential(1)==2

        # do a finite derivative
        drapprox = 1.e-5
        dddudddr = (ddpotential(rcirc+drapprox) - ddpotential(rcirc))/drapprox

    else

        # use the given analytic potential
        dddudddr = dddpotential(rcirc)

    end

    return omega1_circular + (r-rcirc)*((dddudddr + (3/rcirc)*dduddr - (3/rcirc^2)*dudr)/(2*omega1_circular))

end


"""Omega2_expansion(dψ/dr,d²ψ/dr²,r,rcirc)
first-order Taylor expansion for the azimuthal frequency
"""
function Omega2_expansion(dpotential::Function,ddpotential::Function,
                          r::Float64,rcirc::Float64)
    omega2_circular = Omega2_circular(dpotential,rcirc)
    dudr   = dpotential(rcirc)
    dduddr = ddpotential(rcirc)
    h      = r-rcirc
    dO2    = ((1/rcirc)*dduddr - (1/rcirc^2)*(dudr))/(2*omega2_circular)

    return omega2_circular + h*dO2

end


"""Omega1circ_to_radius(Ω₁,dψ/dr,d²ψ/dr²,rmax)
perform backwards mapping from Omega_1 for a circular orbit to radius

can tune [rmin,rmax] for extra optimisation (but not needed)
"""
function Omega1circ_to_radius(omega::Float64,dpotential::Function,ddpotential::Function;rmin::Float64=0.0,rmax::Float64=10000.0)
    r_omega1 = optimize(x -> abs(omega - Omega1_circular(dpotential,ddpotential,x)), rmin    ,rmax  , Brent()).minimizer
    return r_omega1
end


function Omega1circ_to_radius_bisect(omega::Float64,dpotential::Function,ddpotential::Function;rmin::Float64=1.0e-8,rmax::Float64=10000.0,Ziter::Int64=32,verbose::Bool=false)

    # the function to extremise
    extreme(x) = abs(omega - Omega1_circular(dpotential,ddpotential,x))

    # step 1: do a semi-coarse pass over the whole range (24 refinements)
    m,newrmin,newrmax = ExtremiseFunction(extreme,24,rmin,rmax,verbose=verbose,fullreturn=true)

    # step 2: reset the boundaries and do a fine pass
    m = ExtremiseFunction(extreme,Ziter,newrmin,newrmax,verbose=verbose)

    return m

end


"""Omega2circ_to_radius(Ω₂,dψ/dr,d²ψ/dr²[,rmax])
perform backwards mapping from Omega_2 for a circular orbit to radius
"""
function Omega2circ_to_radius(omega::Float64,dpotential::Function,rmax::Float64=1000.)
    r_omega2 = optimize(x -> abs(omega - Omega2_circular(dpotential,x)), 0.    ,rmax  , Brent()).minimizer
    return r_omega2
end


function Omega2circ_to_radius_bisect(omega::Float64,dpotential::Function;rmin::Float64=1.0e-8,rmax::Float64=10000.0,Ziter::Int64=32,verbose::Bool=false)

    # the function to extremise
    extreme(x) = abs(omega - Omega2_circular(dpotential,x))

    # step 1: do a semi-coarse pass over the whole range (24 refinements)
    # the coarse refinement level depends on how small of radii we want to get to
    m,newrmin,newrmax = ExtremiseFunction(extreme,24,rmin,rmax,verbose=verbose,fullreturn=true)

    # step 2: reset the boundaries and do a fine pass
    m = ExtremiseFunction(extreme,Ziter,newrmin,newrmax,verbose=verbose)

    return m

end



"""make_betac(dψ/dr,d²ψ/dr²[,numr,Omega0])
do a high-resolution interpolation to get \beta_c(\alpha), the frequency O2/O1 frequency ratio as a function of O1.

@IMPROVE: find Omega0 adaptively
@IMPROVE: make the interpolation range adaptive in radius
@IMPROVE: decide on best mapping for interplation range (currently log)
@IMPROVE: remove anonymous functions for alpha_c,beta_c

"""
function make_betac(dpotential::Function,ddpotential::Function,numr::Int64=2000,Omega0::Float64=1.)

    # define the circular frequencies
    alpha_c(x) = Omega1_circular(dpotential,ddpotential,x)    # alpha_c(r)
    beta_c(x)  = Omega2_circular(dpotential,x)/alpha_c(x)     # beta_c(r)

    # @IMPROVE: make this range adaptive
    testu = 10 .^ LinRange(5.,-5.,numr)

    # initialise blank arrays
    garr = Array{Float64}(undef, (numr))
    farr = Array{Float64}(undef, (numr))

    # fill the arrays
    for u = 1:numr
        garr[u] = alpha_c(testu[u])/Omega0
        farr[u] = beta_c(testu[u])
    end

    # compute a linear interpolation of the arrays
    beta_c = LinearInterpolation(garr,farr)

    #return beta_c

    # convert this to a function for easier Julia compile
    calculate_betac(x::Float64) = beta_c(x)
    return calculate_betac

end

"""beta_circ(alpha_circ, dψ/dr,d²ψ/dr²[Omega0, rmax])
return \beta_c(\alpha), the frequency O2/O1 frequency ratio as a function of O1.

@IMPROVE: find Omega0 adaptively

"""
function beta_circ(alpha_circ::Float64,dpotential::Function,ddpotential::Function,Omega0::Float64=1.,rmax::Float64=1000.)

    # define the circular frequencies
    omega1 = Omega0 * alpha_circ
    #rcirc = Omega1circ_to_radius(omega1,dpotential,ddpotential,rmax)
    rcirc = Omega1circ_to_radius_bisect(omega1,dpotential,ddpotential,rmax)

    omega2 = Omega2_circular(dpotential,rcirc)

    return omega2/omega1

end
