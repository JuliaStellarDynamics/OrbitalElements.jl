


"""
first-order expansion of energy equation near a circular orbit
"""
function Ecirc_expansion(potential::Function,dpotential::Function,ddpotential::Function,r::Float64,rcirc::Float64)
    # compute the Taylor expansion of Ecirc function at r instead of rcirc
    ur     = potential(rcirc)
    dudr   = dpotential(rcirc)
    dduddr = ddpotential(rcirc)
    h      = r-rcirc
    return (0.5*rcirc*dudr + ur) + h*(0.5*r*dduddr + 1.5*dudr)
end

"""
first-order expansion of angular momentum equation near a circular orbit
"""
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

"""
energy as a function of rp,ra (or a,e)
"""
function E_from_rpra_pot(potential::Function,dpotential::Function,ddpotential::Function,
                         rp::Float64,ra::Float64;
                         TOLECC::Float64=ELTOLECC)


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

"""
angular momentum as a function of rp,ra (or a,e)
"""
function L_from_rpra_pot(potential::Function,dpotential::Function,ddpotential::Function,
                         rp::Float64,ra::Float64;
                         TOLECC::Float64=ELTOLECC)


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

"""
combined energy + angular momentum as a function of rp,ra (or a,e)
"""
function EL_from_rpra_pot(potential::Function,dpotential::Function,ddpotential::Function,
                          rp::Float64,ra::Float64;
                          TOLECC::Float64=ELTOLECC)


    # check the tolerance
    a,ecc = ae_from_rpra(rp,ra)

    if ecc<TOLECC
        # switch to the expanded case
        return Ecirc_expansion(potential,dpotential,ddpotential,rp,a),Lcirc_expansion(potential,dpotential,ddpotential,rp,a)
    else
        # the analytic version of the angular momentum
        # Fouvry 21 eq A5
        return ((ra^2)*potential(ra) - (rp^2)*potential(rp))/(ra^2 - rp^2),sqrt(2*(potential(ra) - potential(rp))/(rp^(-2) - ra^(-2)))
    end

end

"""
combined energy + angular momentum as a function of (a,e)
"""
function EL_from_ae_pot(potential::Function,dpotential::Function,ddpotential::Function,
                        a::Float64,ecc::Float64;
                        TOLECC::Float64=ELTOLECC)



    rp,ra = rpra_from_ae(a,ecc)

    # check the tolerance
    if ecc<TOLECC
        # switch to the expanded case
        return Ecirc_expansion(potential,dpotential,ddpotential,rp,a),Lcirc_expansion(potential,dpotential,ddpotential,rp,a)
    else
        # the analytic version of the angular momentum
        # Fouvry 21 eq A5
        return ((ra^2)*potential(ra) - (rp^2)*potential(rp))/(ra^2 - rp^2),sqrt(2*(potential(ra) - potential(rp))/(rp^(-2) - ra^(-2)))
    end

end


"""
energy and angular momentum derivatives for a given rp,ra
"""
function dEdL_from_rpra_pot(potential::Function,dpotential::Function,ddpotential::Function,
                            rp::Float64,ra::Float64;
                            da::Float64=0.0001,de::Float64=0.0001,
                            TOLECC::Float64=ELTOLECC)


    # check the central tolerance
    a,ecc = ae_from_rpra(rp,ra)

    if ecc<TOLECC
        # switch to the expanded case
        ecc = TOLECC
    end

    # compute the varying values
    rph,rah = rpra_from_ae(a+da,ecc)
    rpr,rar = rpra_from_ae(a,ecc+de)

    # the central values
    Ec = ((ra^2)*potential(ra) - (rp^2)*potential(rp))/(ra^2 - rp^2)
    Lc = sqrt(2*(potential(ra) - potential(rp))/(rp^(-2) - ra^(-2)))

    # the da values
    Eh = ((rah^2)*potential(rah) - (rph^2)*potential(rph))/(rah^2 - rph^2)
    Lh = sqrt(2*(potential(rah) - potential(rph))/(rph^(-2) - rah^(-2)))

    # the de values
    Er = ((rar^2)*potential(rar) - (rpr^2)*potential(rpr))/(rar^2 - rpr^2)
    Lr = sqrt(2*(potential(rar) - potential(rpr))/(rpr^(-2) - rar^(-2)))

    dEda = (Eh-Ec)/(da)
    dEde = (Er-Ec)/(de)
    dLda = (Lh-Lc)/(da)
    dLde = (Lr-Lc)/(de)

    return Ec,Lc,dEda,dEde,dLda,dLde

end
