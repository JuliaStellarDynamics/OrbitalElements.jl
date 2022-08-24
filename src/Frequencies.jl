#=

wrapper to select the methodology for computing the frequencies

only selecting Henon anomaly mapping for now,
but this is where one could select for different anomalies

=#

# bring in the frequency mapping
include("Henon/Frequencies.jl")
# bring in the anomaly mapping (i.e. f(u))
include("Henon/Ufunc.jl")

# bring in the frequency inversion
include("Utils/NumericalInversion.jl")


########################################################################
#
# (a,e) -> (Ω1,Ω2) mapping : Wrappers
#
########################################################################

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e)
"""
function ComputeFrequenciesAE(ψ::Function,
                              dψ::Function,
                              d2ψ::Function,
                              d3ψ::Function,
                              d4ψ::Function,
                              a::Float64,e::Float64;
                              action::Bool=false,
                              TOLECC::Float64=0.001,
                              verbose::Int64=0,
                              NINT::Int64=32,
                              EDGE::Float64=0.01)

    return HenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;action=action,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)
end

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,a,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e)
EXCEPT fourth derivative
"""
function ComputeFrequenciesAE(ψ::Function,
                              dψ::Function,
                              d2ψ::Function,
                              a::Float64,e::Float64;
                              action::Bool=false,
                              TOLECC::Float64=0.001,
                              verbose::Int64=0,
                              NINT::Int64=32,
                              EDGE::Float64=0.01,
                              FDIFF::Float64=1.e-8)

    # define a numerical fourth derivative
    d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/FDIFF

    return ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;action=action,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)
end

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,a,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e)
EXCEPT third derivative
"""
function ComputeFrequenciesAE(ψ::Function,
                              dψ::Function,
                              d2ψ::Function,
                              a::Float64,e::Float64;
                              action::Bool=false,
                              TOLECC::Float64=0.001,
                              verbose::Int64=0,
                              NINT::Int64=32,
                              EDGE::Float64=0.01,
                              FDIFF::Float64=1.e-8)

    # define a numerical third derivative
    d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF
    # Nul fourth derivative
    d4ψ(x::Float64) = 0.

    return ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;action=action,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)
end


########################################################################
#
# (a,e) -> (Ω1,Ω2) mapping : derivatives wrappers
#
########################################################################

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψa,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
"""
function ComputeFrequenciesAEWithDeriv(ψ::Function,
                                       dψ::Function,
                                       d2ψ::Function,
                                       d3ψ::Function,
                                       d4ψ::Function,
                                       a::Float64,
                                       e::Float64;
                                       da::Float64=0.0001,
                                       de::Float64=0.0001,
                                       TOLECC::Float64=0.001,
                                       verbose::Int64=0,
                                       NINT::Int64=32,
                                       EDGE::Float64=0.01)

        # first, check for values that need to be expanded

        # grid is structured like
        # (Ω1h,Ω2h) [+da]
        #    ^
        # (Ω1c,Ω2c)-> (Ω1r,Ω2r) [+de]

        # @IMPROVE watch out for close to TOLECC, will fail across boundary
        Ω1c,Ω2c = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)

        # the offset in a
        Ω1h,Ω2h = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a+da,e;action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)

        # the offset in e
        # if this is already a radial orbit, don't go to super radial
        if ecc+de > 1.0
            de *= -1.0
        end

        Ω1r,Ω2r = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e+de;action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)

        dΩ1da = (Ω1h-Ω1c)/da
        dΩ2da = (Ω2h-Ω2c)/da

        dΩ1de = (Ω1r-Ω1c)/de
        dΩ2de = (Ω2r-Ω2c)/de

        return Ω1c,Ω2c,dΩ1da,dΩ2da,dΩ1de,dΩ2de
end

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψa,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
EXCEPT fourth derivative
"""
function ComputeFrequenciesAEWithDeriv(ψ::Function,
                                       dψ::Function,
                                       d2ψ::Function,
                                       d3ψ::Function,
                                       a::Float64,
                                       e::Float64;
                                       da::Float64=0.0001,
                                       de::Float64=0.0001,
                                       TOLECC::Float64=0.001,
                                       verbose::Int64=0,
                                       NINT::Int64=32,
                                       FDIFF::Float64=1.e-8)

    # define a numerical fourth derivative
    d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/FDIFF

    return ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;da=da,de=de,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)
end

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψa,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
EXCEPT third derivative
"""
function ComputeFrequenciesAEWithDeriv(ψ::Function,
                                       dψ::Function,
                                       d2ψ::Function,
                                       a::Float64,
                                       e::Float64;
                                       da::Float64=0.0001,
                                       de::Float64=0.0001,
                                       TOLECC::Float64=0.001,
                                       verbose::Int64=0,
                                       NINT::Int64=32,
                                       FDIFF::Float64=1.e-8)

    # define a numerical third derivative
    d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF
    # Nul fourth derivative
    d4ψ(x::Float64) = 0.

    return ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;da=da,de=de,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)
end


########################################################################
#
# (Ω1,Ω2) -> (a,e) mapping : Wrappers
#
########################################################################

"""ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,a,ecc[,eps,maxiter,TOLECC,TOLA])
wrapper to select which type of inversion to compute for (Omega1,Omega2)->(a,e)
"""
function ComputeAEFromFrequencies(ψ::Function,
                                  dψ::Function,
                                  d2ψ::Function,
                                  d3ψ::Function,
                                  Ω1::Float64,Ω2::Float64,
                                  eps::Float64=1*10^(-12),
                                  maxiter::Int64=1000,
                                  TOLECC::Float64=0.001,TOLA::Float64=0.0001,
                                  da::Float64=0.0001,de::Float64=0.0001,
                                  verbose::Int64=0)

        # use adaptive da, de branches
        # da max(0.0001,0.01a)
        # de min(max(0.0001,0.1a*e)

        a,e,iter,finaltol = AEFromOmega1Omega2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ;
                                                    eps=eps,ITERMAX=maxiter,
                                                    TOLECC=TOLECC,TOLA=TOLA,da=da,de=de,
                                                    verbose=verbose)

        return a,e
end


########################################################################
#
# (E,L) -> (α,β) mapping : Jacobian
#
########################################################################

"""
compute the jacobian J = |d(E,L)/d(α,β)| = |d(E,L)/d(a,e)|/|d(α,β)/d(a,e)|
"""
function JacELToAlphaBetaAE(ψ::Function,
                            dψ::Function,
                            d2ψ::Function,
                            d3ψ::Function,
                            d4ψ::Function,
                            a::Float64,
                            e::Float64;
                            NINT::Int64=64,
                            EDGE::Float64=0.02,
                            Omega0::Float64=1.0,
                            TOLECC::Float64=0.001)


    # the (E,L) -> (a,e) Jacobian (in Utils/ComputeEL.jl)
    #Jac_EL_AE = JacELToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLECC=TOLECC)
    Jac_EL_AE = JacELToAE(ψ,dψ,d2ψ,a,e,TOLECC=TOLECC)

    # the (alpha,beta) -> (a,e) Jacobian (below)
    Jac_AB_AE = JacAlphaBetaToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT=NINT,EDGE=EDGE,Omega0=Omega0)

    # compute the Jacobian
    Jac = Jac_EL_AE/Jac_AB_AE

    # do some cursory checks for quality
    if Jac < 0.0
        return 0.0
    end

    if isnan(Jac)
        return 0.0
    end

    return Jac

end

"""

@ATTENTION can use the isochrone-specific if you are using an isochrone. Otherwise this is a bit costly.


"""
function JacAlphaBetaToAE(ψ::Function,
                          dψ::Function,
                          d2ψ::Function,
                          d3ψ::Function,
                          d4ψ::Function,
                          a::Float64,
                          e::Float64;
                          NINT::Int64=64,
                          EDGE::Float64=0.02,
                          Omega0::Float64=1.0)

    # calculate the frequency derivatives
    α,β,∂α∂a,∂α∂e,∂β∂a,∂β∂e = OrbitalElements.DHenonThetaFreqRatiosAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT=NINT,EDGE=EDGE,Omega0=Omega0)

    # return the Jacobian
    Jacαβae = abs(∂α∂a*∂β∂e - ∂β∂a*∂α∂e)

end

"""

@ATTENTION this takes (a,e) as arguments.
@ATTENTION this combines several numerical derivatives; please take care!

@IMPROVE add massaging parameters for numerical derivatives
@IMPROVE fix boundary values when using limited development
@IMPROVE noisy at the boundaries
"""
function JacELToAlphaBetaAE(a::Float64,
                            ecc::Float64,
                            ψ::Function,
                            dψ::Function,
                            d2ψ::Function,
                            Ω₀::Float64=1.0;
                            nancheck::Bool=false,
                            NINT::Int64=64)

    tmpecc = ecc
    # to be fixed for limited development...
    if ecc>0.99
        tmpecc=0.99
    end

    if ecc<0.01
        #println("faking the eccentricity...")
        tmpecc=0.01
    end

    # get all numerical derivatives

    # these are dangerous, and break down fairly easily.
    Ω1c,Ω2c,dΩ1da,dΩ2da,dΩ1de,dΩ2de = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,a,tmpecc,NINT=NINT)

    # this is nearly always safe
    # the (E,L) -> (a,e) Jacobian (in Utils/ComputeEL.jl)
    #Jac_EL_AE = JacELToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)
    Jac_EL_AE = JacELToAE(ψ,dψ,d2ψ,a,tmpecc)


    J_o1o2_ae = abs(dΩ1da*dΩ2de - dΩ1de*dΩ2da)

    # check for NaN or zero values
    if nancheck
        if isnan(Jac_EL_AE)
            println("OrbitalElements.Frequencies.JacELToAlphaBetaAE: J_EL_ae is NaN for a=$a,e=$ecc")
            return 0.0
        end

        if Jac_EL_AE <= 0.0
            println("OrbitalElements.Frequencies.JacELToAlphaBetaAE: J_EL_ae is 0 for a=$a,e=$ecc")
            return 0.0
        end

        if isnan(J_o1o2_ae)
            println("OrbitalElements.Frequencies.JacELToAlphaBetaAE: J_o12_ae is NaN for a=$a,e=$ecc")
            return 0.0
        end

        if J_o1o2_ae <= 0.0
            println("OrbitalElements.Frequencies.JacELToAlphaBetaAE: J_o12_ae is 0 for a=$a,e=$ecc")
            return 0.0
        end
    end

    # combine and return
    return Ω1c*Ω₀*Jac_EL_AE/J_o1o2_ae

end