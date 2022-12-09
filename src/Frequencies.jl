#=

wrapper to select the methodology for computing the frequencies

only selecting Henon anomaly mapping for now,
but this is where one could select for different anomalies

=#

# bring in the frequency mapping
include("Henon/Frequencies.jl")
# bring in the anomaly mapping (i.e. f(u))
include("Henon/Ufunc.jl")


########################################################################
#
# (a,e) -> (Ω1,Ω2) mapping : Wrappers
#
########################################################################

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e[,TOLECC,VERBOSE])
wrapper to select which type of frequency computation to perform, from (a,e)
"""
@inline function ComputeFrequenciesAE(ψ::Function,
                              dψ::Function,
                              d2ψ::Function,
                              d3ψ::Function,
                              d4ψ::Function,
                              a::Float64,
                              e::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64}

    return HenonΘFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end

function ComputeFrequenciesJAE(ψ::Function,
                              dψ::Function,
                              d2ψ::Function,
                              d3ψ::Function,
                              d4ψ::Function,
                              a::Float64,
                              e::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64,Float64}

    return HenonΘFrequenciesJAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e[,TOLECC,VERBOSE])
wrapper to select which type of frequency computation to perform, from (a,e)
EXCEPT fourth derivative
"""
function ComputeFrequenciesAE(ψ::Function,
                              dψ::Function,
                              d2ψ::Function,
                              d3ψ::Function,
                              a::Float64,e::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64}

    # define a numerical fourth derivative
    d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/params.FDIFF

    return ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e[,TOLECC,VERBOSE])
wrapper to select which type of frequency computation to perform, from (a,e)
EXCEPT third derivative
"""
function ComputeFrequenciesAE(ψ::Function,
                              dψ::Function,
                              d2ψ::Function,
                              a::Float64,e::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64}

    # define a numerical third derivative
    d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF

    # Nul fourth derivative
    d4ψ(x::Float64) = 0.

    return ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
end


########################################################################
#
# (a,e) -> (Ω1,Ω2) mapping : derivatives wrappers
#
########################################################################

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e[,da,de,TOLECC,VERBOSE,NINT,EDGE])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES

Presently Henon-specific
"""
@inline function ComputeFrequenciesAEWithDeriv(ψ::Function,
                                       dψ::Function,
                                       d2ψ::Function,
                                       d3ψ::Function,
                                       d4ψ::Function,
                                       a::Float64,
                                       e::Float64,
                                       params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

        # first, check for values that need to be expanded

        da, de = params.da, params.de
        # grid is structured like
        # (Ω1h,Ω2h) [+da]
        #    ^
        # (Ω1c,Ω2c)-> (Ω1r,Ω2r) [+de]

        # @IMPROVE watch out for close to TOLECC, will fail across boundary
        #Ω1c,Ω2c = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;action=false,TOLECC=TOLECC,VERBOSE=VERBOSE,NINT=NINT,EDGE=EDGE)
        Ω1c,Ω2c = HenonΘFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

        # the offset in a
        #Ω1h,Ω2h = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a+da,e;action=false,TOLECC=TOLECC,VERBOSE=VERBOSE,NINT=NINT,EDGE=EDGE)
        a2 = a+da
        Ω1h,Ω2h = HenonΘFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a2,e,params)
        # the offset in e
        # if this is already a radial orbit, don't go to super radial
        e2 = e+de
        if e2 > 1.0
            de *= -1.0
            e2 = e+de
        end

        #Ω1r,Ω2r = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e+de;action=false,TOLECC=TOLECC,VERBOSE=VERBOSE,NINT=NINT,EDGE=EDGE)
        Ω1r,Ω2r = HenonΘFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e2,params)

        dΩ1da = (Ω1h-Ω1c)/da
        dΩ2da = (Ω2h-Ω2c)/da

        dΩ1de = (Ω1r-Ω1c)/de
        dΩ2de = (Ω2r-Ω2c)/de

        return Ω1c,Ω2c,dΩ1da,dΩ2da,dΩ1de,dΩ2de
end

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,a,e[,da,de,TOLECC,VERBOSE,NINT,FDIFF,EDGE])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
EXCEPT fourth derivative
"""
function ComputeFrequenciesAEWithDeriv(ψ::Function,
                                       dψ::Function,
                                       d2ψ::Function,
                                       d3ψ::Function,
                                       a::Float64,
                                       e::Float64,
                                       params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

    # define a numerical fourth derivative
    d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/params.FDIFF

    return ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

end

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,a,e[,TOLECC,VERBOSE])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
EXCEPT third derivative
"""
function ComputeFrequenciesAEWithDeriv(ψ::Function,
                                       dψ::Function,
                                       d2ψ::Function,
                                       a::Float64,
                                       e::Float64,
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
@inline function ComputeActionsAE(ψ::Function,
                              dψ::Function,
                              d2ψ::Function,
                              d3ψ::Function,
                              a::Float64,
                              e::Float64,
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
@inline function ComputeActionsAEWithDeriv(ψ::Function,
                                       dψ::Function,
                                       d2ψ::Function,
                                       d3ψ::Function,
                                       a::Float64,
                                       e::Float64,
                                       params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}

        # first, check for values that need to be expanded

        da, de = params.da, params.de
        # grid is structured like
        # (Jh,Lh) [+da]
        #    ^
        # (Jc,Lc)-> (Jr,Lr) [+de]

        # @IMPROVE watch out for close to TOLECC, will fail across boundary
        Jc, Lc = ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,params)

        # the offset in a
        a2 = a+da
        Jh, Lh = ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a2,e,params)
        # the offset in e
        # if this is already a radial orbit, don't go to super radial
        e2 = e+de
        if e2 > 1.0
            de *= -1.0
            e2 = e+de
        end

        Jr,Lr = ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e2,params)

        dJda = (Jh-Jc)/da
        dLda = (Lh-Lc)/da

        dJde = (Jr-Jc)/de
        dLde = (Lr-Lc)/de

        return Jc, Lc, dJda, dLda, dJde, dLde
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
@inline function ComputeAEFromFrequencies(ψ::Function,
                                  dψ::Function,
                                  d2ψ::Function,
                                  d3ψ::Function,
                                  d4ψ::Function,
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
@inline function ComputeAEFromActions(ψ::Function,
                                  dψ::Function,
                                  d2ψ::Function,
                                  d3ψ::Function,
                                  d4ψ::Function,
                                  J::Float64,L::Float64,
                                  params::OrbitsParameters)::Tuple{Float64,Float64}

        a, e, _, _ = AEFromJLBrute(J,L,ψ,dψ,d2ψ,d3ψ,d4ψ,params)

        return a, e
end



########################################################################
#
# (E,L) -> (α,β) mapping : Jacobian
#
########################################################################

"""
compute the jacobian J = |d(E,L)/d(α,β)| = |d(E,L)/d(a,e)|/|d(α,β)/d(a,e)|
"""
function JacELToαβAE(ψ::Function,
                     dψ::Function,
                     d2ψ::Function,
                     d3ψ::Function,
                     d4ψ::Function,
                     a::Float64,
                     e::Float64,
                     params::OrbitsParameters)::Float64


    # the (E,L) -> (a,e) Jacobian (in Utils/ComputeEL.jl)
    Jac_EL_AE = JacELToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    # the (α,β) -> (a,e) Jacobian (below)
    Jac_AB_AE = JacαβToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    # compute the Jacobian
    Jac = Jac_EL_AE/Jac_AB_AE

    # do some cursory checks for quality
    if (Jac < 0.0) || isinf(Jac) || isnan(Jac)
        Jac = 0.0
    end

    return Jac
end

"""

@ATTENTION can use the isochrone-specific if you are using an isochrone. Otherwise this is a bit costly.


"""
function JacαβToAE(ψ::Function,
                          dψ::Function,
                          d2ψ::Function,
                          d3ψ::Function,
                          d4ψ::Function,
                          a::Float64,
                          e::Float64,
                          params::OrbitsParameters)::Float64

    # calculate the frequency derivatives
    _, _, ∂α∂a, ∂α∂e, ∂β∂a, ∂β∂e = OrbitalElements.DHenonΘFreqRatiosAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

    # return the Jacobian
    return abs(∂α∂a*∂β∂e - ∂β∂a*∂α∂e)

end