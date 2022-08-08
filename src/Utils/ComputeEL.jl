"""
energy as a function of (a,e) for a given potential ψ (and its derivatives)
signature with third derivative of the potential
"""
function EFromAE(ψ::Function,
                 dψ::Function,
                 d2ψ::Function,
                 d3ψ::Function,
                 a::Float64,
                 e::Float64;
                 TOLECC::Float64=ELTOLECC)

    if e<TOLECC
        # switch to the expanded case
        return EcircExpansion(ψ,dψ,d2ψ,d3ψ,a,e)
    else
        # the analytic version of the energy
        return ((1+e)^(2)*ψ(a*(1+e)) - (1-e)^(2)*ψ(a*(1-e))) / (4e)
    end
end

"""
energy as a function of (rp,ra) for a given potential ψ (and its derivatives)
signature with third derivative of the potential
"""
function EFromRpRa(ψ::Function,
                   dψ::Function,
                   d2ψ::Function,
                   d3ψ::Function,
                   rp::Float64,
                   ra::Float64;
                   TOLECC::Float64=ELTOLECC)

    a,e = ae_from_rpra(rp,ra)

    return EFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLECC=ELTOLECC)
end

"""
energy as a function of (rp,ra) for a given potential ψ (and its derivatives)
signature with NO third derivative of the potential
"""
function EFromRpRa(ψ::Function,
                   dψ::Function,
                   d2ψ::Function,
                   rp::Float64,
                   ra::Float64;
                   TOLECC::Float64=ELTOLECC,
                   FDIFF::Float64=1.e-8)

    a,e = ae_from_rpra(rp,ra)

    # define a numerical third derivative
    d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF

    return EFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLECC=ELTOLECC)
end


"""
angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)
"""
function LFromAE(ψ::Function,
                 dψ::Function,
                 d2ψ::Function,
                 d3ψ::Function,
                 a::Float64,
                 e::Float64;
                 TOLECC::Float64=ELTOLECC)

    if e<TOLECC
        # switch to the expanded case
        return LcircExpansion(ψ,dψ,d2ψ,d3ψ,a,e)
    else
        # the analytic version of the angular momentum
        return a * (1-(e)^(2)) * sqrt( (ψ(a*(1+e)) - ψ(a*(1-e))) / (2e) )
    end
end

"""
angular momentum as a function of (rp,ra) for a given potential ψ (and its derivatives)
signature with third derivative of the potential
"""
function LFromRpRa(ψ::Function,
                   dψ::Function,
                   d2ψ::Function,
                   d3ψ::Function,
                   rp::Float64,
                   ra::Float64;
                   TOLECC::Float64=ELTOLECC)

    a,e = ae_from_rpra(rp,ra)

    return LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLECC=ELTOLECC)
end

"""
angular momentum as a function of (rp,ra) for a given potential ψ (and its derivatives)
signature with NO third derivative of the potential
"""
function LFromRpRa(ψ::Function,
                   dψ::Function,
                   d2ψ::Function,
                   rp::Float64,
                   ra::Float64;
                   TOLECC::Float64=ELTOLECC,
                   FDIFF::Float64=1.e-8)

    a,e = ae_from_rpra(rp,ra)

    # define a numerical third derivative
    d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF

    return LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLECC=ELTOLECC)
end

"""
combined energy + angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)
"""
function ELFromAE(ψ::Function,
                        dψ::Function,
                        d2ψ::Function,
                        d3ψ::Function,
                        a::Float64,
                        e::Float64;
                        TOLECC::Float64=ELTOLECC)

    E = EFromAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)
    L = LFromAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)

    return E, L
end

"""
combined energy + angular momentum as a function of (rp,ra) for a given potenial ψ (and its derivatives)
"""
function ELFromRpRa(ψ::Function,
                        dψ::Function,
                        d2ψ::Function,
                        d3ψ::Function,
                        a::Float64,
                        e::Float64;
                        TOLECC::Float64=ELTOLECC)

    a,e = ae_from_rpra(rp,ra)
    E = EFromAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)
    L = LFromAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)

    return E, L
end

"""
Second-order expansion of energy equation near a circular orbit
"""
function EcircExpansion(ψ::Function,
                         dψ::Function,
                         d2ψ::Function,
                         d3ψ::Function,
                         a::Float64,
                         e::Float64)

    # compute the Taylor expansion of E
    return (0.5*a*dψ(a) + ψ(a)) + (0.5*a*dψ(a) + 0.5*(a)^(2)*d2ψ(a) + (a)^(3)*d3ψ(a)/12) * (e)^(2)
end

"""
Second-order expansion of angular momentum equation near a circular orbit
"""
function LcircExpansion(ψ::Function,
                         dψ::Function,
                         d2ψ::Function,
                         d3ψ::Function,
                         a::Float64,
                         e::Float64)

    # compute the Taylor expansion of L
    Lcirc = (sqrt(a))^(3)*sqrt(dψ(a))
    return Lcirc + ((a)^(5)*d3ψ(a)/(12*Lcirc) - Lcirc) * (e)^(2)
end

"""
energy and angular momentum derivatives w.r.t. (a,e)
"""
function dELFromAE(ψ::Function,
                   dψ::Function,
                   d2ψ::Function,
                   d3ψ::Function,
                   d4ψ::Function,
                   a::Float64,
                   e::Float64;
                   TOLECC::Float64=ELTOLECC)

    E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)

    if e<TOLECC
        # switch to the expanded case
        ∂E∂a, ∂E∂e, ∂L∂a, ∂L∂e = dELcircExpansion(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)
        return E, L, ∂E∂a, ∂E∂e, ∂L∂a, ∂L∂e

    else

        # the analytic version of the energy and angular momentum derivatives w.r.t. (a,e)
        rp, ra = rpra_from_ae(a,e)
        ψrp, ψra, dψrp, dψra = ψ(rp), ψ(ra), dψ(rp), dψ(ra)

        # Difference between potential at apocenter and pericenter
        ψdiff = ψra - ψrp

        ∂E∂a = ((1+e)^(3)*dψra - (1-e)^(3)*dψrp) / (4e)
        ∂E∂e = (((e)^(2)-1)*ψdiff + a*e*(1+e)^(2)*dψra + a*e*(1-e)^(2)*dψrp) / (4*(e)^(2))

        dLdenom = 2*sqrt(2*e*ψdiff)

        ∂L∂a = (1-(e)^(2)) * (2*ψdiff + ra*dψra - rp*dψrp)  / (dLdenom)
        ∂L∂e = - a * ((1+3*(e)^(2))*ψdiff - a*e*(1-(e)^(2))*(dψra + dψrp)) / (e*dLdenom)

        return E, L, ∂E∂a, ∂E∂e, ∂L∂a, ∂L∂e
    end
end

"""
Second-order of energy and angular momentum derivatives w.r.t. (a,e) near circular orbits.
"""
function dELcircExpansion(ψ::Function,
                          dψ::Function,
                          d2ψ::Function,
                          d3ψ::Function,
                          d4ψ::Function,
                          a::Float64,
                          e::Float64)

    # compute all potential and derivative values
    ψa, dψa, ddψa, dddψa, ddddψa = ψ(a), dψ(a), d2ψ(a), d3ψ(a), d4ψ(a)

    # compute taylor expansion of partial derivatives for E
    ∂E∂a = 0.5*(3*dψa + a*ddψa) + (0.5*dψa + 1.5*a*ddψa + 0.75*(a)^(2)*dddψa + (a)^(3)*ddddψa/12) * (e)^(2)
    ∂E∂e = (a*dψa + (a)^(2)*ddψa + (a)^(3)*dddψa/6) * e

    sqa   = sqrt(a)
    sqdψa = sqrt(dψa)

    # compute taylor expansion of partial derivatives for L
    ∂L∂a = 0.5*(3*sqrt(a)*sqdψa + (sqa)^(3)*ddψa/sqdψa) - 0.5 * sqa * (3*(dψa)^(2) + a*dψa*ddψa - 7*(a)^(2)*dψa*dddψa/12 + (a)^(3)*(ddψa*dddψa - 2*dψa*ddddψa)/12) / ((sqdψa)^(3)) * (e)^(2)
    ∂L∂e = (sqa)^(3) * ((a)^(2)*dddψa/(6*sqdψa) - 2*sqdψa) * e

    return ∂E∂a, ∂E∂e, ∂L∂a, ∂L∂e
end
