

"""
    αβFromUV(u,v,n1,n2,ωmin,ωmax)

mapping from (α,β) to (u,v)

Fouvry & Prunet B5

This signature specifies ωmin and ωmax, to avoid extra calls.

"""
function αβFromUV(u::Float64,v::Float64,
                  n1::Int64,n2::Int64,
                  ωmin::Float64,ωmax::Float64)::Tuple{Float64,Float64}

    if n2 == 0
        β  = v
        α = (1.0/(2*n1))*((ωmax-ωmin)*u + ωmin + ωmax)
    else
        α = v
        β  = (1.0/(n2*v))*(0.5*((ωmax-ωmin)*u + ωmin + ωmax) - n1*v)
    end

    return α,β
end


"""
    UVFromαβ(α,β,n1,n2,ωmin,ωmax)

mapping from  (u,v) to (α,β), from Fouvry & Prunet Appendix B

@IMPROVE, this has rounding error: concern?

OrbitalElements.UVFromαβ(0.5,0.7,-3,4,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr)

"""
function UVFromαβ(α::Float64,β::Float64,
                  n1::Int64,n2::Int64,
                  ωmin::Float64,ωmax::Float64)::Tuple{Float64,Float64}

    # Equation B1
    ωval = n1*α + n2*β*α

    # Equation B3
    u = (2*ωval - ωmax - ωmin)/(ωmax-ωmin)

    # Equation B4
    v = (n2==0) ? β : α

    return u,v
end

"""
    RenormalizedJacUVToαβ(n1,n2,u,v)

Renormalized jacobian of the (u,v) ↦ (α,β) mapping, i.e. 2/(ωmax-ωmin) * |∂(α,β)/∂(u,v)|

@ATTENTION: depends on the resonance number (as the mapping does).
"""
function RenormalizedJacUVToαβ(n1::Int64,n2::Int64,
                               u::Float64,v::Float64)::Float64
                               
    return (n2==0) ? 1.0 / abs(n1) : 1.0 / abs(n2*v)
end
