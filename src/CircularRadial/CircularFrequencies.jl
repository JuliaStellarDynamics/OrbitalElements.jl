
using Optim                  # for zero-finding

# for circular orbits, we can use exact relations with potential derivatives
function Omega1_circular(dpotential::Function,ddpotential::Function,r::Float64)
    return sqrt(ddpotential(r) + 3*dpotential(r)/r)
end


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


function Omega2_circular(dpotential::Function,r::Float64)
    return sqrt(dpotential(r)/r)
end

function Omega2_expansion(dpotential::Function,ddpotential::Function,
                          r::Float64,rcirc::Float64)
    omega2_circular = Omega2_circular(dpotential,rcirc)
    dudr   = dpotential(rcirc)
    dduddr = ddpotential(rcirc)
    h      = r-rcirc
    dO2    = ((1/rcirc)*dduddr - (1/rcirc^2)*(dudr))/(2*omega2_circular)

    return omega2_circular + h*dO2

end


function Omega1circ_to_radius(omega::Float64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.)
    # perform the backwards mapping from Omega_1 for a circular orbit to radius
    r_omega1 = optimize(x -> abs(omega - Omega1_circular(dpotential,ddpotential,x)), 0.    ,rmax  , Brent()).minimizer
    return r_omega1
end

function Omega2circ_to_radius(omega::Float64,dpotential::Function,rmax::Float64=1000.)
    # perform the backwards mapping from Omega_2 for a circular orbit to radius
    r_omega2 = optimize(x -> abs(omega - Omega2_circular(dpotential,x)), 0.    ,rmax  , Brent()).minimizer
    return r_omega2
end
