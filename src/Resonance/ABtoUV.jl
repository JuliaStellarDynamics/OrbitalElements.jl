"""AlphaBetaFromUV(u,v,n1,n2,dψ,d2ψ,rmax,Omega0)

mapping from (alpha,beta) to (u,v)

Fouvry & Prunet B5

@IMPROVE, this has rounding error: concern?
@IMPROVE, Omega0 isnt really optional, so we should perhaps not set a default?

"""
function AlphaBetaFromUV(u::Float64,v::Float64,
                         n1::Int64,n2::Int64,
                         dψ::Function,d2ψ::Function,
                         rmax::Float64=1000.,Omega0::Float64=1.)


    ωmin,ωmax = FindWminWmax(n1,n2,dψ,d2ψ,rmax,Omega0)

    return AlphaBetaFromUV(u,v,n1,n2,ωmin,ωmax)

end

"""AlphaBetaFromUV(u,v,n1,n2,ωmin,ωmax)

mapping from (alpha,beta) to (u,v)

Fouvry & Prunet B5

This signature specifies ωmin and ωmax, to avoid extra calls.

"""
function AlphaBetaFromUV(u::Float64,v::Float64,
                         n1::Int64,n2::Int64,
                         ωmin::Float64,ωmax::Float64)

    if n2 == 0
        beta  = v
        alpha = (1.0/(2*n1))*((ωmax-ωmin)*u + ωmin + ωmax)
    else
        alpha = v
        beta  = (1.0/(n2*v))*(0.5*((ωmax-ωmin)*u + ωmin + ωmax) - n1*v)
    end

    return alpha,beta
end

"""UVFromAlphaBeta(alpha,beta,n1,n2,dψ,d2ψ[,rmax,Omega0])

mapping from  (u,v) to (alpha,beta)

@IMPROVE, this has rounding error: concern?

OrbitalElements.UVFromAlphaBeta(0.5,0.7,-3,4,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr)

"""
function UVFromAlphaBeta(alpha::Float64,beta::Float64,
                           n1::Int64,n2::Int64,
                           dψ::Function,d2ψ::Function,
                           rmax::Float64=1000.,Omega0=1.)

    ωmin,ωmax = FindWminWmax(n1,n2,dψ,d2ψ,rmax,Omega0)

    return UVFromAlphaBeta(alpha,beta,n1,n2,ωmin,ωmax)

end

"""UVFromAlphaBeta(alpha,beta,n1,n2,ωmin,ωmax )

mapping from  (u,v) to (alpha,beta), from Fouvry & Prunet Appendix B

@IMPROVE, this has rounding error: concern?

OrbitalElements.UVFromAlphaBeta(0.5,0.7,-3,4,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr)

"""
function UVFromAlphaBeta(alpha::Float64,beta::Float64,
                           n1::Int64,n2::Int64,
                           ωmin::Float64,ωmax::Float64)

    # Equation B1
    wval = n1*alpha + n2*beta*alpha

    # Equation B3
    u = (2*wval - ωmax - ωmin)/(ωmax-ωmin)

    # Equation B4
    if (n2==0)
        v = beta
    else
        v = alpha
    end

    return u,v

end

"""
using the definitions for (alpha, beta) and (u,v), compute the Jacobian.
@ATTENTION, to match eq. B6, this has the 2/(ωmax-ωmin) term already absorbed into it. therefore, not formally the Jacobian, but adds the dimensional removal.
"""
function JacalphabetaToUV(n1::Int64,n2::Int64,w_min::Float64,w_max::Float64,v::Float64)

    if n2 ==0
        return (2.0/(w_max-w_min))*abs((w_max-w_min) * (1/(2*n1)))
    else
        return (2.0/(w_max-w_min))*abs((w_max-w_min) * (1/(2*n2*v)))
    end
end
