

"""αβFromUV(u,v,n1,n2,ωmin,ωmax)

mapping from (α,β) to (u,v)

Fouvry & Prunet B5

This signature specifies ωmin and ωmax, to avoid extra calls.

"""
function αβFromUV(u::Float64,v::Float64,
                  n1::Int64,n2::Int64,
                  ωmin::Float64,ωmax::Float64)

    if n2 == 0
        β  = v
        α = (1.0/(2*n1))*((ωmax-ωmin)*u + ωmin + ωmax)
    else
        α = v
        β  = (1.0/(n2*v))*(0.5*((ωmax-ωmin)*u + ωmin + ωmax) - n1*v)
    end

    return α,β
end


"""UVFromαβ(α,β,n1,n2,ωmin,ωmax )

mapping from  (u,v) to (α,β), from Fouvry & Prunet Appendix B

@IMPROVE, this has rounding error: concern?

OrbitalElements.UVFromαβ(0.5,0.7,-3,4,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr)

"""
function UVFromαβ(α::Float64,β::Float64,
                  n1::Int64,n2::Int64,
                  ωmin::Float64,ωmax::Float64)

    # Equation B1
    ωval = n1*α + n2*β*α

    # Equation B3
    u = (2*ωval - ωmax - ωmin)/(ωmax-ωmin)

    # Equation B4
    if (n2==0)
        v = β
    else
        v = α
    end

    return u,v
end

"""
using the definitions for (α, β) and (u,v), compute the Jacobian.
@ATTENTION, to match eq. B6, this has the 2/(ωmax-ωmin) term already absorbed into it. therefore, not formally the Jacobian, but adds the dimensional removal.
"""
function JacαβToUV(n1::Int64,n2::Int64,ωmin::Float64,ωmax::Float64,v::Float64)

    if (n2==0)
        return (2.0/(ωmax-ωmin))*abs((ωmax-ωmin) * (1/(2n1)))
    else
        return (2.0/(ωmax-ωmin))*abs((ωmax-ωmin) * (1/(2n2*v)))
    end
end
