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



"""ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e)
"""
function ComputeFrequenciesAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                              a::Float64,ecc::Float64;
                              action::Bool=false,
                              TOLECC::Float64=0.001,
                              verbose::Int64=0,
                              NINT::Int64=32,
                              EDGE::Float64=0.01)

    out = HenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,ecc,action=true,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)

    if action
        # out is a tuple with three values: O1,O2,Jr
        Ω1,Ω2,Jr = out
        return Ω1,Ω2,Jr
    else
        # out is a tuple with two values: O1,O2
        Ω1,Ω2 = out
        return Ω1,Ω2
    end

end



"""ComputeFrequenciesAE(ψ,dψ,d2ψ,a,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e)
"""
function ComputeFrequenciesAE(ψ::Function,dψ::Function,d2ψ::Function,
                              a::Float64,ecc::Float64;
                              action::Bool=false,
                              TOLECC::Float64=0.001,
                              verbose::Int64=0,
                              NINT::Int64=32)

    out = ComputeFrequenciesHenonAE(ψ,dψ,d2ψ,a,ecc,action=true,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

    if action
        Ω1,Ω2,a1 = out
        return Ω1,Ω2,a1
    else
        Ω1,Ω2 = out
        return Ω1,Ω2
    end

end


"""compute_ae_from_frequencies(ψ,dψ,d2ψ,a,ecc[,eps,maxiter,TOLECC,TOLA])
wrapper to select which type of inversion to compute for (Omega1,Omega2)->(a,e)
"""
function compute_ae_from_frequencies(ψ::Function,dψ::Function,d2ψ::Function,
                                     Ω1::Float64,Ω2::Float64,
                                     eps::Float64=1*10^(-12),
                                     maxiter::Int64=1000,
                                     TOLECC::Float64=0.001,TOLA::Float64=0.0001,
                                     da::Float64=0.0001,de::Float64=0.0001,
                                     verbose::Int64=0)

        # use adaptive da, de branches
        # da max(0.0001,0.01a)
        # de min(max(0.0001,0.1a*e)

        a,e,iter,finaltol = ae_from_omega1omega2_brute(Ω1,Ω2,ψ,dψ,d2ψ,eps,maxiter,TOLECC,TOLA,da,de,verbose)

        ntries = 0
        while (iter == maxiter+1) | (iter <= 0)
            # double the da step to scan through space
            da = 2da
            a,e,iter,finaltol = ae_from_omega1omega2_brute(Ω1,Ω2,ψ,dψ,d2ψ,eps,maxiter,TOLECC,TOLA,da,de,verbose)
            ntries += 1
            if ntries > 3
                break
            end
        end

        # more optional massages to try and go smaller

        return a,e
end


"""ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,a,ecc[,eps,maxiter,TOLECC,TOLA])
wrapper to select which type of inversion to compute for (Omega1,Omega2)->(a,e)
"""
function ComputeAEFromFrequencies(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                                  Ω1::Float64,Ω2::Float64,
                                  eps::Float64=1*10^(-12),
                                  maxiter::Int64=1000,
                                  TOLECC::Float64=0.001,TOLA::Float64=0.0001,
                                  da::Float64=0.0001,de::Float64=0.0001,
                                  verbose::Int64=0)

        # use adaptive da, de branches
        # da max(0.0001,0.01a)
        # de min(max(0.0001,0.1a*e)

        a,e,iter,finaltol = AEFromOmega1Omega2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d2ψ,eps,maxiter,TOLECC,TOLA,da,de,verbose)

        return a,e
end


"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψa,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES

@IMPROVE: could add action derivatives here? more copacetic with analytic derivatives anyway
"""
function ComputeFrequenciesRpRaWithDeriv(ψ::Function,
                                         dψ::Function,
                                         d2ψ::Function,
                                         rp::Float64,
                                         ra::Float64;
                                         dr::Float64=1.e-6,
                                         TOLECC::Float64=0.001,
                                         verbose::Int64=0,
                                         NINT::Int64=32)


        d3ψ(r::Float64)::Float64 = (d2ψ(r+dr) - d2ψ(r))/dr
        # grid is structured like
        # (Ω1h,Ω2h) [+da]
        #    ^
        # (Ω1c,Ω2c)-> (Ω1r,Ω2r) [+de]

        # @IMPROVE watch out for close to TOLECC, will fail across boundary
        a,e = ae_from_rpra(rp,ra)
        Ω1c,Ω2c = HenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e,action=false,TOLECC=TOLECC,NINT=NINT)


        # need to check 'polarity': that is, if a circular orbit, don't do rp+dr; do rp-dr
        if (rp+dr > ra)
            dr *= -1.0
        end

        #ah,eh = 0.5*(rp+dr+ra),(ra-(rp+dr))/(rp+ra+dr)
        ah,eh = ae_from_rpra(rp+dr,ra)
        Ω1h,Ω2h = HenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,ah,eh,action=false,TOLECC=TOLECC,NINT=NINT)
        dΩ1drp = (Ω1h-Ω1c)/dr
        dΩ2drp = (Ω2h-Ω2c)/dr

        # need to check 'polarity': that is, if a circular orbit, don't do ra-dr; do ra+dr
        if (ra+dr < rp)
            dr *= -1.0
        end

        ar,er = ae_from_rpra(rp,ra+dr)
        Ω1r,Ω2r = HenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,ar,er,action=false,TOLECC=TOLECC,NINT=NINT)
        dΩ1dra = (Ω1r-Ω1c)/dr
        dΩ2dra = (Ω2r-Ω2c)/dr

        return Ω1c,Ω2c,dΩ1drp,dΩ2drp,dΩ1dra,dΩ2dra
end


"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψa,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
INCLUDING third potential derivative
"""
function ComputeFrequenciesAEWithDeriv(ψ::Function,
                                       dψ::Function,
                                       d2ψ::Function,
                                       d3ψ::Function,
                                       a::Float64,
                                       ecc::Float64;
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
        Ω1c,Ω2c = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,ecc,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)

        # the offset in a
        Ω1h,Ω2h = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a+da,ecc,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)

        # the offset in e
        # if this is already a radial orbit, don't go to super radial
        if ecc+de > 1.0
            de *= -1.0
        end

        Ω1r,Ω2r = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,ecc+de,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)

        dΩ1da = (Ω1h-Ω1c)/da
        dΩ2da = (Ω2h-Ω2c)/da

        dΩ1de = (Ω1r-Ω1c)/de
        dΩ2de = (Ω2r-Ω2c)/de

        return Ω1c,Ω2c,dΩ1da,dΩ2da,dΩ1de,dΩ2de
end



"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψa,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES

@IMPROVE: could add action derivatives here? more copacetic with analytic derivatives anyway
"""
function ComputeFrequenciesAEWithDeriv(ψ::Function,
                                       dψ::Function,
                                       d2ψ::Function,
                                       a::Float64,
                                       ecc::Float64;
                                       da::Float64=0.0001,
                                       de::Float64=0.0001,
                                       TOLECC::Float64=0.001,
                                       verbose::Int64=0,
                                       NINT::Int64=32)

        # first, check for values that need to be expanded

        # grid is structured like
        # (Ω1h,Ω2h) [+da]
        #    ^
        # (Ω1c,Ω2c)-> (Ω1r,Ω2r) [+de]

        # @IMPROVE watch out for close to TOLECC, will fail across boundary
        Ω1c,Ω2c = ComputeFrequenciesHenonAE(ψ,dψ,d2ψ,a,ecc,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        Ω1h,Ω2h = ComputeFrequenciesHenonAE(ψ,dψ,d2ψ,a+da,ecc,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        # if this is already a radial orbit, don't go to super radial
        if ecc+de > 1.0
            de *= -1.0
        end

        Ω1r,Ω2r = ComputeFrequenciesHenonAE(ψ,dψ,d2ψ,a,ecc+de,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        dΩ1da = (Ω1h-Ω1c)/da
        dΩ2da = (Ω2h-Ω2c)/da

        dΩ1de = (Ω1r-Ω1c)/de
        dΩ2de = (Ω2r-Ω2c)/de

        return Ω1c,Ω2c,dΩ1da,dΩ2da,dΩ1de,dΩ2de
end


function ComputeFrequenciesAEWithDerivCircular(ψ::Function,
                                               dψ::Function,
                                               d2ψ::Function,
                                               a::Float64,
                                               ecc::Float64,
                                               da::Float64=0.0001,
                                               de::Float64=0.001,
                                               TOLECC::Float64=0.001,
                                               verbose::Int64=0,
                                               NINT::Int64=32)



        # grid is structured like
        # (Ω1h,Ω2h) [+da]
        #    ^
        # (Ω1c,Ω2c)-> (Ω1r,Ω2r) [+de]

        # get the frequencies using the epicyclic approximation
        Ω1c = Omega1_circular(dψ,d2ψ,a)
        Ω2c = Omega2_circular(dψ,a)

        # use the epicyclic approximation to step forward a tiny bit in radius
        Ω1h = Omega1_circular(dψ,d2ψ,a+da)
        Ω2h = Omega2_circular(dψ,a+da)

        # take a fairly large step in eccentricity to make sure we reach a safe zone
        # @ATTENTION, TOLECC needs to be smaller than de for this to work
        Ω1r,Ω2r = ComputeFrequenciesHenonAE(ψ,dψ,d2ψ,a,de,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        dΩ1da = (Ω1h-Ω1c)/da
        dΩ2da = (Ω2h-Ω2c)/da

        dΩ1de = (Ω1r-Ω1c)/de
        dΩ2de = (Ω2r-Ω2c)/de

        return Ω1c,Ω2c,dΩ1da,dΩ2da,dΩ1de,dΩ2de
end



"""
compute the jacobian J = |d(E,L)/d(alpha,beta)| = |d(E,L)/d(a,e)|/|d(alpha,beta)/d(a,e)|
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



"""ComputeFrequenciesRpRa(ψ,dψ,d2ψ,r_peri,r_apo[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e)
"""
function ComputeFrequenciesRpRa(ψ::Function,dψ::Function,d2ψ::Function,
                                  r_peri::Float64,r_apo::Float64,TOLECC::Float64=0.001,verbose::Int64=0,NINT=32)

        Ω1,Ω2 = ComputeFrequenciesHenonRpRa(ψ,dψ,d2ψ,r_peri,r_apo,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        return Ω1,Ω2
end
