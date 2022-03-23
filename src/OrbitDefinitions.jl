#
# generic orbit transformations

function ae_from_rpra(rp::Float64,ra::Float64)
    #=
    function to translate pericentre and apocentre to semi-major axis and eccentricity
    =#
    return (rp+ra)/2,(ra-rp)/(rp+ra)
end

function rpra_from_ae(a::Float64,e::Float64)
    #=
    function to translate semi-major axis and eccentricity to pericentre and apocentre
    =#
    return a*(1-e),a*(1+e)
end



function find_j(potential::Function,r::Float64,kappa::Float64)
    # compute the angular momentum for a given orbit at radius r, with given kappa
    dudr = dpotential_numerical(potential,r)
    jmax = sqrt(r*r*r*dudr);
    J = jmax*kappa;
    return J
end

function find_K(potential::Function,r::Float64,J::Float64)
    # compute kappa for a given orbit at radius r, with given angular momentum
    jmax = Lcirc_numerical(potential::Function,r::Float64)
    kappa = J/jmax;
    # would like to do some sort of expansion here
    #if kappa>1.0
    #    kappa=2-kappa
    #end
    if kappa<0.0
        kappa=0.0
    end
    return kappa
end


function Ecirc(potential::Function,r::Float64,E::Float64)
    # compute the energy of a circular orbit at some radius
    # must define the potential and potential derivative a priori
    ur = potential(r)
    dudr = dpotential_numerical(potential,r)
    return  abs(E - 0.5*r*dudr - ur)
end


function denom(potential::Function,r::Float64,E::Float64,J::Float64)
    # the main function to root-find, this is the effective potential
    ur = potential(r)
    return abs(2.0*(E-ur)*r*r - J*J)
end


function make_orbit(potential::Function,E::Float64,K::Float64,rmax::Float64=100000.)
    # initialise an orbit in E,K space
    # will not check that E is valid for the model!
    # rmax must be defined in order to recover appropriate roots
    r_circ = optimize(r -> Ecirc(potential,r,E), 0.    ,rmax  , Brent()).minimizer
    J      = find_j(potential,r_circ,K)
    r_apo  = optimize(r -> denom(potential,r,E,J), r_circ,rmax  , Brent()).minimizer
    r_peri = optimize(r -> denom(potential,r,E,J), 0.    ,r_circ, Brent()).minimizer
    return r_peri,r_apo,r_circ,J
end

function make_orbit_ae(potential::Function,a::Float64,ecc::Float64,rmax::Float64=100000.,TOLECC::Float64=0.00005)
    # initialise an orbit in a,e space
    r_peri,r_apo = rpra_from_ae(a,ecc)
    E = E_from_rpra_pot(potential,r_peri,r_apo,TOLECC)
    J = L_from_rpra_pot(potential,r_peri,r_apo,TOLECC)
    # will not check that E is valid for the model!
    # rmax must be defined in order to recover appropriate roots
    if ecc<TOLECC
        r_circ = a#0.5*(r_peri+r_apo)
    else
        r_circ = optimize(r -> Ecirc(potential,r,E), 0.    ,rmax  , Brent()).minimizer
    end
    return r_peri,r_apo,r_circ,J
end
