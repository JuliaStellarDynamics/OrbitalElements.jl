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


"""MakeOrbitAE

initialise an orbit in (a,e) space, returning rp, ra, rcirc, L

@IMPROVE: will not check that E is valid for the model!
@IMPROVE: something better than defaulting to circular below some eccentricity?

"""
function MakeOrbitAE(ψ::Function,dψ::Function,d2ψ::Function,
                       a::Float64,e::Float64;
                       TOLECC::Float64=0.00005)

    rp, ra = RpRafromAE(a,e)

    # get (E,L)
    E = EFromRpRa(ψ,dψ,d2ψ,rp,ra;TOLECC=TOLECC)
    L = LFromRpRa(ψ,dψ,d2ψ,rp,ra;TOLECC=TOLECC)

    if ecc<TOLECC
        rc = a
    else
        rc = bisection(r -> E - Ecirc(ψ,dψ,r),rp,ra)
    end

    return rp, ra, rc, L

end
