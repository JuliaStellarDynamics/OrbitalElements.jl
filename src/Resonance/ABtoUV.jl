"""alphabeta_from_uv

mapping from (alpha,beta) to (u,v)

Fouvry & Prunet B5

@IMPROVE, this has rounding error: concern?

"""
function alphabeta_from_uv(u::Float64,v::Float64,
                           n1::Int64,n2::Int64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.,Omega0::Float64=1.)


    wmin,wmax = find_wmin_wmax(n1,n2,dpotential,ddpotential,rmax,Omega0)

    if n2 == 0
        beta  = v
        alpha = (1/(2n1))*((wmax-wmin)*u + wmin + wmax)
    else
        alpha = v
        beta  = (1/(n2*v))*(0.5*((wmax-wmin)*u + wmin + wmax) - n1*v)
    end

    return alpha,beta
end

"""uv_from_alphabeta(alpha,beta,n1,n2,dpotential,ddpotential[,rmax,Omega0])

mapping from  (u,v) to (alpha,beta)

@IMPROVE, this has rounding error: concern?

OrbitalElements.uv_from_alphabeta(0.5,0.7,-3,4,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr)

"""
function uv_from_alphabeta(alpha::Float64,beta::Float64,
                           n1::Int64,n2::Int64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.,Omega0=1.)

    wmin,wmax = find_wmin_wmax(n1,n2,dpotential,ddpotential,rmax,Omega0)

    wval = n1*alpha + n2*beta*alpha

    u = (2*wval - wmax - wmin)/(wmax-wmin)

    if (n2==0)
        v = beta
    else
        v = alpha
    end

    return u,v

end

"""
using the definitions for (alpha, beta) and (u,v), compute the Jacobian.
"""
function Jacalphabeta_to_uv(n1::Int64,n2::Int64,w_min::Float64,w_max::Float64,v::Float64)

    if n2 ==0
        return abs(w_max-w_min) * (1/(2n1))
    else
        return abs(w_max-w_min) * (1/(2n2*v))
    end
end
