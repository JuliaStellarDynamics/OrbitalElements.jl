"""
basic orbit transformations

"""

"""AEfromRpRa

function to translate pericentre and apocentre to semi-major axis and eccentricity

"""
@inline function AEfromRpRa(rp::Float64,ra::Float64)::Tuple{Float64,Float64}

    return (rp+ra)/2,(ra-rp)/(rp+ra)
end

"""RpRafromAE

function to translate semi-major axis and eccentricity to pericentre and apocentre

"""
@inline function RpRafromAE(a::Float64,e::Float64)::Tuple{Float64,Float64}

    return a*(1-e),a*(1+e)
end



"""Ecirc

compute the energy of a circular orbit at some radius
must define the potential and potential derivative a priori

"""
@inline function Ecirc(ψ::Function,dψ::Function,r::Float64)::Float64

    return  ψ(r) + 0.5*r*dψ(r)
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
@inline function Vrad(ψ::Function,
              dψ::Function,
              d2ψ::Function,
              d3ψ::Function,
              u::Float64,
              a::Float64,
              e::Float64,
              params::OrbitsParameters)::Float64

    E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)

    r = ru(u,a,e)

    vrSQ = 2*(E - ψeff(ψ,r,L))

    if vrSQ < 0.0
        return 0.0
    else
        return sqrt(vrSQ)
    end
end
