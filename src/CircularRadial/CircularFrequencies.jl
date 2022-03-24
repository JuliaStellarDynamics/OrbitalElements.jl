#=
CircularFrequencies.jl

Special treatment for circular orbits
  for circular orbits, we can use exact relations with potential derivatives

=#

"""Omega1_circular

radial frequency for circular orbits, from the epicyclic approximation

inputs
-------------
- dpotential   : (Function) the potential
- ddpotential  : (Function) the potential derivative
- r            : (Float64)  the radius at which to evaluate

returns
-------------
- freq1        : (Float64) the radial frequency from the epicyclic approximation

"""
function Omega1_circular(dpotential::Function,ddpotential::Function,r::Float64)

    return sqrt(ddpotential(r) + 3*dpotential(r)/r)
end

"""
Omega2_circular

azimuthal frequency for circular orbits, from the epicyclic approximation

inputs
-------------
- ddpotential  : (Function) the potential derivative
- r            : (Float64)  the radius at which to evaluate

returns
-------------
- freq2        : (Float64) the azimuthal frequency from the epicyclic approximation

"""
function Omega2_circular(dpotential::Function,r::Float64)

    return sqrt(dpotential(r)/r)
end

"""
Omega1_expansion

first-order Taylor expansion for the radial frequency


"""
function Omega1_expansion(dpotential::Function,ddpotential::Function,
                          r::Float64,rcirc::Float64,
                          dddpotential::Function=f(x)=x+1)

    omega1_circular = Omega1_circular(dpotential,ddpotential,r)
    dudr = dpotential(rcirc)
    dduddr = ddpotential(rcirc)
    h    = r-rcirc
    if dddpotential(1)==2
        dddudddr = 0.
    else
        dddudddr = dddpotential(rcirc)
    end

    return omega1_circular + (r-rcirc)*((dddudddr + (3/rcirc)*dduddr - (3/rcirc^2)*dudr)/(2*omega1_circular))

end


"""
Omega2_expansion

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


"""
Omega1circ_to_radius

perform backwards mapping from Omega_1 for a circular orbit to radius

"""
function Omega1circ_to_radius(omega::Float64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.)
    r_omega1 = optimize(x -> abs(omega - Omega1_circular(dpotential,ddpotential,x)), 0.    ,rmax  , Brent()).minimizer
    return r_omega1
end

"""
Omega2circ_to_radius

perform backwards mapping from Omega_2 for a circular orbit to radius

"""
function Omega2circ_to_radius(omega::Float64,dpotential::Function,rmax::Float64=1000.)
    r_omega2 = optimize(x -> abs(omega - Omega2_circular(dpotential,x)), 0.    ,rmax  , Brent()).minimizer
    return r_omega2
end


"""make_betac

do a high-resolution interpolation to get \beta_c(alpha), the frequency O2/O1 frequency ratio as a function of O1.

"""
function make_betac(dpotential::Function,ddpotential::Function,numr::Int64=2000,Omega0::Float64=1.)


    alpha_c(x) = Omega1_circular(dpotential,ddpotential,x)    # alpha_c(r)
    beta_c(x)  = Omega2_circular(dpotential,x)/alpha_c(x)     # beta_c(r)


    # @IMPROVE: make this range adaptive
    testu = 10 .^ LinRange(5.,-5.,numr)

    garr = Array{Float64}(undef, (numr))
    farr = Array{Float64}(undef, (numr))

    for u = 1:numr
        garr[u] = alpha_c(testu[u])/Omega0
        farr[u] = beta_c(testu[u])
    end

    beta_c = LinearInterpolation(garr,farr)

    return beta_c
end
