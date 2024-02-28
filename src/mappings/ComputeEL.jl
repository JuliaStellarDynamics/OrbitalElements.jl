"""Definitions for handling energy and angular momentum

Strategies:
-compute energy and angular momentum from the definitions for (a,e)
-compute derivatives of energy and angular momentum
-add switches for near-circular orbits (E,L,dE,dL)
-include Jacobian to transform between (E,L) and (a,e)
-auxiliary functions to do same calculations as a function of (rp,ra)

"""



########################################################################
#
# (a,e) -> (E,L) mapping : generic case
#
########################################################################

"""
    EFromAE(ψ,dψ,a,e,params)

energy as a function of (a,e) for a given potential ψ (and its derivatives)
"""
function EFromAE(model::Potential,
                 a::Float64,e::Float64,
                 params::OrbitalParameters=OrbitalParameters())::Float64

    # Handling edges interpolations
    # IMPORTANT : has to be first !
    fun(atemp::Float64,etemp::Float64) = EFromAE(model,atemp,etemp,params)
    E = EdgeHandle(fun,a,e,params)
    if !(isnothing(E))
        return E
    end

    # Edge cases
    if (a == 0.)
        return ψ(model,0.)
    elseif (e == 0.)
        return Ecirc(model,a)
    elseif (e == 1.)
        return Erad(model,a)
    end

    # Generic
    return ((1+e)^(2)*ψ(model,a*(1+e)) - (1-e)^(2)*ψ(model,a*(1-e))) / (4e)
end

"""
    LFromAE(ψ,dψ,a,e,params)

angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)
"""
function LFromAE(model::Potential,
                 a::Float64,e::Float64,
                 params::OrbitalParameters=OrbitalParameters())::Float64

    # Handling edges interpolations
    # IMPORTANT : has to be first !
    fun(atemp::Float64,etemp::Float64) = LFromAE(model,atemp,etemp,params)
    L = EdgeHandle(fun,a,e,params)
    if !(isnothing(L))
        return L
    end

    # Edge cases
    if (a == 0.)
        return 0.
    elseif (e == 0.)
        return Lcirc(model,a)
    elseif (e == 1.)
        return 0.
    end

    # Generic
    return a * (1-(e)^(2)) * sqrt( (ψ(model,a*(1+e)) - ψ(model,a*(1-e))) / (2e) )
end

"""
    ELFromAE(ψ,dψ,a,e,params)

combined energy + angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)
"""
function ELFromAE(model::Potential,
                  a::Float64,e::Float64,
                  params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64}

    E = EFromAE(model,a,e,params)
    L = LFromAE(model,a,e,params)

    return E, L
end

########################################################################
#
# (a,e) -> (E,L) mapping : radial values
#
########################################################################

"""
    Erad(ψ,a)

energy for an exactly radial orbit
"""
function Erad(model::Potential,
              a::Float64)::Float64

    return ψ(model,2*a)
end

########################################################################
#
# (a,e) -> (E,L) mapping : circular values
#
########################################################################

"""
    Ecirc(ψ,dψ,a)

energy for an exactly circular orbit.
"""
function Ecirc(model::Potential,
               a::Float64)::Float64

    return ψ(model,a) + 0.5*a*dψ(model,a)
end

"""
    Lcirc(dψ,a)

angular momentum for an exactly circular orbit.
"""
function Lcirc(model::Potential,
               a::Float64)::Float64

    return (sqrt(a))^(3)*sqrt(dψ(model,a))
end


########################################################################
#
# (a,e) -> (E,L) mapping : derivatives, generic case
#
########################################################################

"""
    ComputeELAEWithDeriv(ψ,dψ,a,e,params)
"""
function ComputeELAEWithDeriv(model::Potential,
                              a::Float64,e::Float64,
                              params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

    # Function to differentiate
    fun(atemp::Float64,etemp::Float64) = ELFromAE(model,atemp,etemp,params)
    # Perform differentiation
    floc, ∂f∂a, ∂f∂e = NumericalDerivativeAE(fun,a,e,params)
    # Recast results
    E, L = floc
    ∂E∂a, ∂L∂a = ∂f∂a
    ∂E∂e, ∂L∂e = ∂f∂e

    return E, L, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e
end


########################################################################
#
# (a,e) -> (E,L) mapping : Jacobian
#
########################################################################

"""
    JacAEToEL(ψ,dψ,d2ψ,a,e,params)

Jacobian of the (a,e) ↦ (E,L) mapping, i.e. |∂(E,L)/∂(a,e)|
"""
function JacAEToEL(model::Potential,
                   a::Float64,e::Float64,
                   params::OrbitalParameters=OrbitalParameters())

    _, _, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e = ComputeELAEWithDeriv(model,a,e,params)

    return abs(∂E∂a*∂L∂e - ∂L∂a*∂E∂e)
end

########################################################################
#
# (E,L) -> (a,e) mapping
#
########################################################################

"""
    ComputeAEFromEL(ψ,dψ,E,L,params)
"""
function ComputeAEFromEL(model::Potential,
                         E::Float64,L::Float64,
                         params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64}

    a, e, _, _ = AEFromELBrute(E,L,model,params)

    return a, e
end


########################################################################
#
# (rp,ra) -> (E,L) mapping : from (a,e) -> (E,L) mapping
#
########################################################################

"""
    EFromRpRa(ψ,dψ,rp,ra,params)

energy as a function of (rp,ra) for a given potential ψ (and its derivatives)
"""
function EFromRpRa(model::Potential,
                   rp::Float64,ra::Float64,
                   params::OrbitalParameters=OrbitalParameters())::Float64

    a,e = AEFromRpRa(rp,ra)

    return EFromAE(model,a,e,params)
end


"""
    LFromRpRa(ψ,dψ,rp,ra,params)

angular momentum as a function of (rp,ra) for a given potential ψ (and its derivatives)
"""
function LFromRpRa(model::Potential,
                   rp::Float64,ra::Float64,
                   params::OrbitalParameters=OrbitalParameters())::Float64

    a,e = AEFromRpRa(rp,ra)

    return LFromAE(model,a,e,params)
end


"""
    ELFromRpRa(ψ,dψ,rp,ra,params)

combined energy + angular momentum as a function of (rp,ra) for a given potenial ψ (and its derivatives)
"""
function ELFromRpRa(model::Potential,
                    rp::Float64,ra::Float64,
                    params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64}

    a,e = AEFromRpRa(rp,ra)

    return ELFromAE(model,a,e,params)
end