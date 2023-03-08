"""
basic orbit transformations

"""

"""AEfromRpRa

function to translate pericentre and apocentre to semi-major axis and eccentricity

"""
function AEFromRpRa(rp::Float64,ra::Float64)::Tuple{Float64,Float64}

    return (rp+ra)/2, (ra-rp)/(rp+ra)
end

"""RpRafromAE

function to translate semi-major axis and eccentricity to pericentre and apocentre

"""
function RpRaFromAE(a::Float64,e::Float64)::Tuple{Float64,Float64}

    return a*(1-e), a*(1+e)
end

"""JacRpRaToAE

Jacobian of the (rp,ra) ↦ (a,e) mapping, i.e. |∂(a,e)/∂(rp,ra)|
"""
function JacRpRaToAE(rp::Float64,ra::Float64)::Float64

    return 1.0 / (rp + ra)
end


"""JacAEToRpRa

Jacobian of the (a,e) ↦ (rp,ra) mapping, i.e. |∂(rp,ra)/∂(a,e)|
"""
function JacAEToRpRa(a::Float64,e::Float64)::Float64

    return 2a
end

########################################################################
#
# vr(u) (action integrand)
#
########################################################################

"""Vrad(ψ,dψ,d2ψd,3ψ,u,a,e[,TOLECC,fun])
vr, radial velocity for computing action
as a function of (a,e)
"""
function Vrad(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
              u::Float64,a::Float64,e::Float64,
              params::OrbitsParameters)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)

    r = ru(u,a,e)

    vrSQ = 2*(E - ψeff(ψ,r,L))

    if (vrSQ < 0.0) || isnan(vrSQ) || isinf(vrSQ)
        return 0.0
    else
        return sqrt(vrSQ)
    end
end
