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
function EFromAE(ψ::F0,dψ::F1,
                 a::Float64,e::Float64,
                 params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function}

    # Handling edges interpolations
    # IMPORTANT : has to be first !
    fun(atemp::Float64,etemp::Float64) = EFromAE(ψ,dψ,atemp,etemp,params)
    E = EdgeHandle(fun,a,e,params)
    if !(isnothing(E))
        return E
    end

    # Edge cases
    if (a == 0.)
        return ψ(0.)
    elseif (e == 0.)
        return Ecirc(ψ,dψ,a)
    elseif (e == 1.)
        return Erad(ψ,a)
    end

    # Generic
    return ((1+e)^(2)*ψ(a*(1+e)) - (1-e)^(2)*ψ(a*(1-e))) / (4e)
end

"""
    LFromAE(ψ,dψ,a,e,params)

angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)
"""
function LFromAE(ψ::F0,dψ::F1,
                 a::Float64,e::Float64,
                 params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function}

    # Handling edges interpolations
    # IMPORTANT : has to be first !
    fun(atemp::Float64,etemp::Float64) = LFromAE(ψ,dψ,atemp,etemp,params)
    L = EdgeHandle(fun,a,e,params)
    if !(isnothing(L))
        return L
    end

    # Edge cases
    if (a == 0.)
        return 0.
    elseif (e == 0.)
        return Lcirc(dψ,a)
    elseif (e == 1.)
        return 0.
    end

    # Generic
    return a * (1-(e)^(2)) * sqrt( (ψ(a*(1+e)) - ψ(a*(1-e))) / (2e) )
end

"""
    ELFromAE(ψ,dψ,a,e,params)

combined energy + angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)
"""
function ELFromAE(ψ::F0,dψ::F1,
                  a::Float64,e::Float64,
                  params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function}

    E = EFromAE(ψ,dψ,a,e,params)
    L = LFromAE(ψ,dψ,a,e,params)

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
function Erad(ψ::Function,
              a::Float64)::Float64

    return ψ(2*a)
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
function Ecirc(ψ::F0,dψ::F1,
               a::Float64)::Float64 where {F0 <: Function, F1 <: Function}

    return ψ(a) + 0.5*a*dψ(a)
end

"""
    Lcirc(dψ,a)

angular momentum for an exactly circular orbit.
"""
function Lcirc(dψ::Function,
               a::Float64)::Float64

    return (sqrt(a))^(3)*sqrt(dψ(a))
end


########################################################################
#
# (a,e) -> (E,L) mapping : derivatives, generic case
#
########################################################################

"""
    ComputeELAEWithDeriv(ψ,dψ,a,e,params)
"""
function ComputeELAEWithDeriv(ψ::F0,dψ::F1,
                              a::Float64,e::Float64,
                              params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function}

    # Function to differentiate
    fun(atemp::Float64,etemp::Float64) = ELFromAE(ψ,dψ,atemp,etemp,params)
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
function JacAEToEL(ψ::F0,dψ::F1,
                   a::Float64,e::Float64,
                   params::OrbitalParameters=OrbitalParameters()) where {F0 <: Function, F1 <: Function}

    _, _, ∂E∂a, ∂L∂a, ∂E∂e, ∂L∂e = ComputeELAEWithDeriv(ψ,dψ,a,e,params)

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
function ComputeAEFromEL(ψ::F0,dψ::F1,
                         E::Float64,L::Float64,
                         params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function}

    a, e, _, _ = AEFromELBrute(E,L,ψ,dψ,params)

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
function EFromRpRa(ψ::F0,dψ::F1,
                   rp::Float64,ra::Float64,
                   params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function}

    a,e = AEFromRpRa(rp,ra)

    return EFromAE(ψ,dψ,a,e,params)
end


"""
    LFromRpRa(ψ,dψ,rp,ra,params)

angular momentum as a function of (rp,ra) for a given potential ψ (and its derivatives)
"""
function LFromRpRa(ψ::F0,dψ::F1,
                   rp::Float64,ra::Float64,
                   params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function}

    a,e = AEFromRpRa(rp,ra)

    return LFromAE(ψ,dψ,a,e,params)
end


"""
    ELFromRpRa(ψ,dψ,rp,ra,params)

combined energy + angular momentum as a function of (rp,ra) for a given potenial ψ (and its derivatives)
"""
function ELFromRpRa(ψ::F0,dψ::F1,
                    rp::Float64,ra::Float64,
                    params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function}

    a,e = AEFromRpRa(rp,ra)

    return ELFromAE(ψ,dψ,a,e,params)
end