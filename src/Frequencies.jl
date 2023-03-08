#=

wrapper to select the methodology for computing the frequencies

only selecting Henon anomaly mapping for now,
but this is where one could select for different anomalies

=#

# bring in the anomaly mapping (i.e. f(u))
include("Henon/Ufunc.jl")
# bring in the frequency mapping
include("Henon/Frequencies.jl")

########################################################################
#
# (α,β) -> (Ω1,Ω2) mapping
#
########################################################################
function αβFromFrequencies(Ω1::Float64,Ω2::Float64,
                           Ω₀::Float64)::Tuple{Float64,Float64}

    return Ω1/Ω₀, Ω2/Ω1
end

function FrequenciesFromαβ(α::Float64,β::Float64,
                           Ω₀::Float64)::Tuple{Float64,Float64}

    return Ω₀*α, Ω₀*α*β
end

function FrequenciesDerivsFromαβDerivs(α::Float64,β::Float64,
                                       dα::Float64,dβ::Float64,
                                       Ω₀::Float64)::Tuple{Float64,Float64}

    return Ω₀*dα, Ω₀*(dα*β + α*dβ) 
end

########################################################################
#
# (a,e) -> (Ω1,Ω2) mapping : Wrappers
#
########################################################################
"""
    αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

wrapper of henon frequency ratio computation.
"""
function αβFromAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                             a::Float64,e::Float64,
                             params::OrbitsParameters)::Tuple{Float64,Float64}

    return αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end

"""
    ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

wrapper to select which type of frequency computation to perform, from (a,e)
"""
function ComputeFrequenciesAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                              a::Float64,e::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64}

    α, β = αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
    return FrequenciesFromαβ(α,β,params.Ω₀)
end

function ComputeFrequenciesJAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                               a::Float64,e::Float64,
                               params::OrbitsParameters)::Tuple{Float64,Float64,Float64}

    Ω1, Ω2 = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
    J = HenonJFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)

    return Ω1, Ω2, J
end

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
wrapper to select which type of frequency computation to perform, from (a,e)
EXCEPT fourth derivative
"""
function ComputeFrequenciesAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                              a::Float64,e::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64}

    # define a numerical fourth derivative
    d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/params.FDIFF

    return ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e,params)
wrapper to select which type of frequency computation to perform, from (a,e)
EXCEPT third derivative
"""
function ComputeFrequenciesAE(ψ::Function,dψ::Function,d2ψ::Function,
                              a::Float64,e::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64}

    # define a numerical third derivative
    d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/params.FDIFF

    # Nul fourth derivative
    d4ψ(x::Float64) = 0.

    return ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end


########################################################################
#
# (a,e) -> (Ω1,Ω2), (α,β) mapping : derivatives wrappers
#
########################################################################

function ComputeαβWithDerivAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              a::Float64,e::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # Numerical derivative points
    ap, da, ep, de = NumDerivPoints(a,e,params.da,params.de,params.TOLA,params.TOLECC)

    # Derivation outside the integral
    α, β = αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    # For a derivatives
    αap, βap = αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,ap,e,params)
    
    # For e derivatives
    αep, βep = αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,ep,params)

    ∂α∂a = (αap-α)/da
    ∂β∂a = (βap-β)/da

    ∂α∂e = (αep-α)/de
    ∂β∂e = (βep-β)/de

    return α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e
end

"""
    ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e[,da,de,TOLECC,VERBOSE,NINT,EDGE])

wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
"""
function ComputeFrequenciesAEWithDeriv(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                                       a::Float64,e::Float64,
                                       params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    Ω₀ = params.Ω₀
    α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    Ω1, Ω2          = FrequenciesFromαβ(α,β,Ω₀)
    ∂Ω1∂a, ∂Ω2∂a    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂a,∂β∂a,Ω₀)
    ∂Ω1∂e, ∂Ω2∂e    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂e,∂β∂e,Ω₀)

    return Ω1, Ω2, ∂Ω1∂a, ∂Ω2∂a, ∂Ω1∂e, ∂Ω2∂e
end

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,a,e[,da,de,TOLECC,VERBOSE,NINT,FDIFF,EDGE])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
EXCEPT fourth derivative
"""
function ComputeFrequenciesAEWithDeriv(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                                       a::Float64,e::Float64,
                                       params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

    # define a numerical fourth derivative
    d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/params.FDIFF

    return ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,a,e[,TOLECC,VERBOSE])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
EXCEPT third derivative
"""
function ComputeFrequenciesAEWithDeriv(ψ::Function,dψ::Function,d2ψ::Function,
                                       a::Float64,e::Float64,
                                       params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

    # define a numerical third derivative
    d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/params.FDIFF
    # Nul fourth derivative
    d4ψ(x::Float64) = 0.

    return ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end

########################################################################
#
# (a,e) -> (J,L) mapping : Wrappers
#
########################################################################

"""ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e[,TOLECC,VERBOSE])
wrapper to select which type of actions computation to perform, from (a,e)
"""
function ComputeActionsAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                          a::Float64,e::Float64,
                          params::OrbitsParameters)::Tuple{Float64,Float64}

    J = HenonJFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
    L = LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
    return J, L
end

########################################################################
#
# (a,e) -> (J,L) mapping : derivatives wrappers
#
########################################################################

"""ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
wrapper to select which type of actions computation to perform, from (a,e), but DERIVATIVES

Presently Henon-specific
"""
function ComputeActionsAEWithDeriv(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                                   a::Float64,e::Float64,
                                   params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # Numerical derivative points
    ap, da, ep, de = NumDerivPoints(a,e,params.da,params.de,params.TOLA,params.TOLECC)

    # Derivation outside the integral
    J, L = ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,params)

    # For a derivatives
    Jap, Lap = ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,ap,e,params)
    
    # For e derivatives
    Jep, Lep = ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,ep,params)

    ∂J∂a = (Jap-J)/da
    ∂L∂a = (Lap-L)/da

    ∂J∂e = (Jep-J)/de
    ∂L∂e = (Lep-L)/de

    return J, L, ∂J∂a, ∂L∂a, ∂J∂e, ∂L∂e
end


########################################################################
#
# (Ω1,Ω2) -> (a,e) mapping : Wrappers
#
########################################################################

# bring in the frequency inversion
include("Utils/NumericalInversion.jl")

"""ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,a,e[,eps,maxiter,TOLECC,TOLA])
wrapper to select which type of inversion to compute for (Omega1,Omega2)->(a,e)
"""
function ComputeAEFromFrequencies(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                                  Ω1::Float64,Ω2::Float64,
                                  params::OrbitsParameters)::Tuple{Float64,Float64}

        # use adaptive da, de branches
        # da max(0.0001,0.01a)
        # de min(max(0.0001,0.1a*e)

        a, e, _, _ = AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ,d4ψ,params)

        return a, e
end


"""ComputeAEFromActions(ψ,dψ,d2ψ,d3ψ,a,e,params)
wrapper to select which type of inversion to compute for (Omega1,Omega2)->(a,e)
"""
function ComputeAEFromActions(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                              J::Float64,L::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64}

        a, e, _, _ = AEFromJLBrute(J,L,ψ,dψ,d2ψ,d3ψ,params)

        return a, e
end



########################################################################
#
# (E,L) -> (α,β) mapping : Jacobian
#
########################################################################

"""
Jacobian of the (α,β) ↦ (E,L) mapping, i.e. |∂(E,L)/∂(α,β)| = |∂(E,L)/∂(a,e)| / |∂(α,β)/∂(a,e)|
"""
function JacαβToELAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                     a::Float64,e::Float64,
                     params::OrbitsParameters)::Float64


    # the (a,e) -> (E,L) Jacobian (in Utils/ComputeEL.jl)
    Jac_AE_To_EL = JacAEToEL(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    # the (α,β) -> (a,e) Jacobian (below)
    Jac_AE_To_αβ = JacAEToαβ(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    # compute the Jacobian
    Jac_αβ_To_EL = Jac_AE_To_EL/Jac_AE_To_αβ

    # do some cursory checks for quality
    if (Jac_αβ_To_EL < 0.0) || isinf(Jac_αβ_To_EL) || isnan(Jac_αβ_To_EL)
        Jac_αβ_To_EL = 0.0
    end

    return Jac_αβ_To_EL
end

"""
Jacobian of the (a,e) ↦ (α,β) mapping, i.e. |∂(α,β)/∂(a,e)|
@ATTENTION can use the isochrone-specific if you are using an isochrone. Otherwise this is a bit costly.
"""
function JacAEToαβ(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                   a::Float64,e::Float64,
                   params::OrbitsParameters)::Float64

    # calculate the frequency derivatives
    _, _, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    # return the Jacobian
    return abs(∂α∂a*∂β∂e - ∂β∂a*∂α∂e)
end