"""
basic orbit transformations

"""

"""AEfromRpRa

function to translate pericentre and apocentre to semi-major axis and eccentricity

"""
function AEfromRpRa(rp::Float64,ra::Float64)

    return (rp+ra)/2,(ra-rp)/(rp+ra)
end

"""RpRafromAE

function to translate semi-major axis and eccentricity to pericentre and apocentre

"""
function RpRafromAE(a::Float64,e::Float64)

    return a*(1-e),a*(1+e)
end


"""Ecirc

compute the energy of a circular orbit at some radius
must define the potential and potential derivative a priori

"""
function Ecirc(ψ::Function,dψ::Function,r::Float64)

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
function Vrad(ψ::Function,
              dψ::Function,
              d2ψ::Function,
              d3ψ::Function,
              u::Float64,
              a::Float64,
              e::Float64;
              TOLECC::Float64=ELTOLECC,
              fun::Function=henon_f)

    E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)

    r = ru(u,a,e;fun=fun)

    vrSQ = 2*(E - ψeff(ψ,r,L))

    if vrSQ < 0.0
        return 0.0
    else
        return sqrt(vrSQ)
    end
end
