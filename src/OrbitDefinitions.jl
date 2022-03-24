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



function Ecirc(potential::Function,dpotential::Function,r::Float64,E::Float64)
    #=
     compute the energy of a circular orbit at some radius
     must define the potential and potential derivative a priori

    =#
    ur   = potential(r)
    dudr = dpotential(potential,r)
    return  abs(E - 0.5*r*dudr - ur)

end


function effective_potential(potential::Function,r::Float64,E::Float64,J::Float64)
    #=

    the main function to root-find, this is the effective potential

    =#
    ur = potential(r)
    return abs(2.0*(E-ur)*r*r - J*J)
end



function make_orbit_ae(potential::Function,a::Float64,ecc::Float64,rmax::Float64=100000.,TOLECC::Float64=0.00005)
    #=

    initialise an orbit in a,e space

    =#

    r_peri,r_apo = rpra_from_ae(a,ecc)

    # get (E,J)
    E = E_from_rpra_pot(potential,r_peri,r_apo,TOLECC)
    J = L_from_rpra_pot(potential,r_peri,r_apo,TOLECC)


    # @IMPROVE: will not check that E is valid for the model!
    # @IMPROVE: rmax must be defined in order to recover appropriate roots
    # @IMPROVE: something better than defaulting to circular below some eccentricity?
    if ecc<TOLECC
        r_circ = a
    else
        r_circ = optimize(r -> Ecirc(potential,dpotential,r,E), 0.    ,rmax  , Brent()).minimizer
    end

    return r_peri,r_apo,r_circ,J

end
