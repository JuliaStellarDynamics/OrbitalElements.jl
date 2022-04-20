"""
basic orbit transformations

"""

"""ae_from_rpra

function to translate pericentre and apocentre to semi-major axis and eccentricity

"""
function ae_from_rpra(rp::Float64,ra::Float64)
    return (rp+ra)/2,(ra-rp)/(rp+ra)
end

"""rpra_from_ae

function to translate semi-major axis and eccentricity to pericentre and apocentre

"""
function rpra_from_ae(a::Float64,e::Float64)
    return a*(1-e),a*(1+e)
end


"""Ecirc

compute the energy of a circular orbit at some radius
must define the potential and potential derivative a priori

"""
function Ecirc(potential::Function,dpotential::Function,r::Float64,E::Float64)
    ur   = potential(r)
    dudr = dpotential(potential,r)
    return  abs(E - 0.5*r*dudr - ur)
end

"""effective_potential

the main function to root-find, this is the effective potential

"""
function effective_potential(potential::Function,r::Float64,E::Float64,L::Float64)

    ur = potential(r)
    return abs(2.0*(E-ur)*r*r - L*L)

end


"""make_orbit_ae

initialise an orbit in (a,e) space, returning rperi,rapo,rcirc,L

@IMPROVE: will not check that E is valid for the model!
@IMPROVE: rmax must be defined in order to recover appropriate roots
@IMPROVE: something better than defaulting to circular below some eccentricity?

"""
function make_orbit_ae(potential::Function,
                       a::Float64,
                       ecc::Float64,
                       rmax::Float64=100000.,
                       TOLECC::Float64=0.00005)

    r_peri,r_apo = rpra_from_ae(a,ecc)

    # get (E,L)
    E = E_from_rpra_pot(potential,r_peri,r_apo,TOLECC)
    L = L_from_rpra_pot(potential,r_peri,r_apo,TOLECC)

    if ecc<TOLECC
        r_circ = a
    else
        r_circ = optimize(r -> Ecirc(potential,dpotential,r,E), 0.    ,rmax  , Brent()).minimizer
    end

    return r_peri,r_apo,r_circ,L

end
