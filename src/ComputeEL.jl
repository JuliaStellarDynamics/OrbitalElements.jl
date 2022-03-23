function Ecirc_expansion(potential::Function,dpotential::Function,ddpotential::Function,r::Float64,rcirc::Float64)
    # compute the Taylor expansion of Ecirc function at r instead of rcirc
    ur     = potential(rcirc)
    dudr   = dpotential(rcirc)
    dduddr = ddpotential(rcirc)
    h      = r-rcirc
    return (0.5*rcirc*dudr + ur) + h*(0.5*r*dduddr + 1.5*dudr)
end

function Lcirc_expansion(potential::Function,dpotential::Function,ddpotential::Function,
                         r::Float64,rcirc::Float64)
    # compute the Taylor expansion of Lcirc function at r instead of rcirc
    ur     = potential(rcirc)
    dudr   = dpotential(rcirc)
    dduddr = ddpotential(rcirc)
    h      = r-rcirc
    dLcirc = ((rcirc^3)*dduddr + 3(rcirc^2)*dudr)/(2*sqrt(rcirc^3*dudr))
    return sqrt((rcirc^3)*dudr) + h*(dLcirc)
end

function E_from_rpra_pot(potential::Function,dpotential::Function,ddpotential::Function,
                         rp::Float64,ra::Float64,
                         TOLECC::Float64=0.005)


    # check the tolerance
    a,ecc = ae_from_rpra(rp,ra)

    if ecc<TOLECC
        # switch to the expanded case
        return Ecirc_expansion(potential,dpotential,ddpotential,rp,a)
    else
        # the analytic version of the energy
        # Fouvry 21 eq A5
        return ((ra^2)*potential(ra) - (rp^2)*potential(rp))/(ra^2 - rp^2)
    end
end

function L_from_rpra_pot(potential::Function,dpotential::Function,ddpotential::Function,
                         rp::Float64,ra::Float64,
                         TOLECC::Float64=0.005)


    # check the tolerance
    a,ecc = ae_from_rpra(rp,ra)

    if ecc<TOLECC
        # switch to the expanded case
        return Lcirc_expansion(potential,dpotential,ddpotential,rp,a)
    else
        # the analytic version of the angular momentum
        # Fouvry 21 eq A5
        return sqrt(2*(potential(ra) - potential(rp))/(rp^(-2) - ra^(-2)))
    end

end
