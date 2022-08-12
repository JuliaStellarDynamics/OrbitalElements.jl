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
                            nancheck::Bool=false)

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
    f1c,f2c,df1da,df2da,df1de,df2de = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,a,tmpecc)

    # this is nearly always safe
    # the (E,L) -> (a,e) Jacobian (in Utils/ComputeEL.jl)
    #Jac_EL_AE = JacELToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)
    Jac_EL_AE = JacELToAE(ψ,dψ,d2ψ,a,tmpecc)


    J_o1o2_ae = abs(df1da*df2de - df1de*df2da)

    # check for NaN or zero values
    if nancheck
        if isnan(Jac_EL_AE)
            println("OrbitalElements.Frequencies.JacELToAlphaBetaAE: J_EL_ae is NaN for a=$a,e=$ecc")
            return 1.0
        end

        if Jac_EL_AE == 0.0
            println("OrbitalElements.Frequencies.JacELToAlphaBetaAE: J_EL_ae is 0 for a=$a,e=$ecc")
            return 1.0
        end

        if isnan(J_o1o2_ae)
            println("OrbitalElements.Frequencies.JacELToAlphaBetaAE: J_o12_ae is NaN for a=$a,e=$ecc")
            return 1.0
        end

        if J_o1o2_ae == 0.0
            println("OrbitalElements.Frequencies.JacELToAlphaBetaAE: J_o12_ae is 0 for a=$a,e=$ecc")
            return 1.0
        end
    end

    # combine and return
    return f1c*Ω₀*Jac_EL_AE/J_o1o2_ae

end



"""compute_frequencies_ae(ψ,dψ/dr,d²ψ/dr²,a,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e)
"""
function compute_frequencies_ae(potential::Function,dpotential::Function,ddpotential::Function,
                                a::Float64,ecc::Float64;action::Bool=false,TOLECC::Float64=0.001,verbose::Int64=0,NINT::Int64=32)

    if action
        f1,f2,a1 = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,ecc,action=true,TOLECC=TOLECC,verbose=verbose,NINT=NINT)
        return f1,f2,a1
    else
        f1,f2 = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,ecc,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)
        return f1,f2
    end
end


"""compute_ae_from_frequencies(ψ,dψ/dr,d²ψ/dr²,a,ecc[,eps,maxiter,TOLECC,TOLA])
wrapper to select which type of inversion to compute for (Omega1,Omega2)->(a,e)
"""
function compute_ae_from_frequencies(potential::Function,dpotential::Function,ddpotential::Function,
                                     omega1::Float64,omega2::Float64,
                                     eps::Float64=1*10^(-12),
                                     maxiter::Int64=1000,
                                     TOLECC::Float64=0.001,TOLA::Float64=0.0001,
                                     da::Float64=0.0001,de::Float64=0.0001,
                                     verbose::Int64=0)

        # use adaptive da, de branches
        # da max(0.0001,0.01a)
        # de min(max(0.0001,0.1a*e)

        a,e,iter,finaltol = ae_from_omega1omega2_brute(omega1,omega2,potential,dpotential,ddpotential,eps,maxiter,TOLECC,TOLA,da,de,verbose)

        ntries = 0
        while (iter == maxiter+1) | (iter <= 0)
            # double the da step to scan through space
            da = 2da
            a,e,iter,finaltol = ae_from_omega1omega2_brute(omega1,omega2,potential,dpotential,ddpotential,eps,maxiter,TOLECC,TOLA,da,de,verbose)
            ntries += 1
            if ntries > 3
                break
            end
        end

        # more optional massages to try and go smaller

        return a,e
end


"""ComputeFrequenciesAEWithDeriv(ψ,dψ/dr,d²ψ/dr²a,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES

@IMPROVE: could add action derivatives here? more copacetic with analytic derivatives anyway
"""
function ComputeFrequenciesRpRaWithDeriv(potential::Function,
                                         dpotential::Function,
                                         ddpotential::Function,
                                         rp::Float64,
                                         ra::Float64;
                                         dr::Float64=0.0001,
                                         TOLECC::Float64=0.001,
                                         verbose::Int64=0,
                                         NINT::Int64=32)



        # grid is structured like
        # (f1h,f2h) [+da]
        #    ^
        # (f1c,f2c)-> (f1r,f2r) [+de]

        # @IMPROVE watch out for close to TOLECC, will fail across boundary
        a,e = 0.5*(rp+ra),(ra-rp)/(rp+ra)
        #f1c,f2c = compute_frequencies_henon(potential,dpotential,ddpotential,rp,ra,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)
        f1c,f2c = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,e,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)
        #println("HARD CHECK $f1c,$f2c,$a,$e")

        # need to check 'polarity': that is, if a circular orbit, don't do rp+dr; do rp-dr
        if (rp+dr > ra)
            dr *= -1.0
        end

        f1h,f2h = compute_frequencies_henon(potential,dpotential,ddpotential,rp+dr,ra,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        df1drp = (f1h-f1c)/dr
        df2drp = (f2h-f2c)/dr

        # need to check 'polarity': that is, if a circular orbit, don't do ra-dr; do ra+dr
        if (ra+dr < rp)
            dr *= -1.0
        end

        f1r,f2r = compute_frequencies_henon(potential,dpotential,ddpotential,rp,ra+dr,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        df1dra = (f1r-f1c)/dr
        df2dra = (f2r-f2c)/dr

        return f1c,f2c,df1drp,df2drp,df1dra,df2dra
end


"""ComputeFrequenciesAEWithDeriv(ψ,dψ/dr,d²ψ/dr²a,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES

@IMPROVE: could add action derivatives here? more copacetic with analytic derivatives anyway
"""
function ComputeFrequenciesAEWithDeriv(potential::Function,
                                       dpotential::Function,
                                       ddpotential::Function,
                                       a::Float64,
                                       ecc::Float64,
                                       da::Float64=0.0001,
                                       de::Float64=0.0001,
                                       TOLECC::Float64=0.001,
                                       verbose::Int64=0,
                                       NINT::Int64=32)

        # first, check for values that need to be expanded

        # grid is structured like
        # (f1h,f2h) [+da]
        #    ^
        # (f1c,f2c)-> (f1r,f2r) [+de]

        # @IMPROVE watch out for close to TOLECC, will fail across boundary
        f1c,f2c = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,ecc,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        f1h,f2h = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a+da,ecc,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        # if this is already a radial orbit, don't go to super radial
        if ecc+de > 1.0
            de *= -1.0
        end

        f1r,f2r = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,ecc+de,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        df1da = (f1h-f1c)/da
        df2da = (f2h-f2c)/da

        df1de = (f1r-f1c)/de
        df2de = (f2r-f2c)/de

        return f1c,f2c,df1da,df2da,df1de,df2de
end


function ComputeFrequenciesAEWithDerivCircular(potential::Function,
                                               dpotential::Function,
                                               ddpotential::Function,
                                               a::Float64,
                                               ecc::Float64,
                                               da::Float64=0.0001,
                                               de::Float64=0.001,
                                               TOLECC::Float64=0.001,
                                               verbose::Int64=0,
                                               NINT::Int64=32)



        # grid is structured like
        # (f1h,f2h) [+da]
        #    ^
        # (f1c,f2c)-> (f1r,f2r) [+de]

        # get the frequencies using the epicyclic approximation
        f1c = Omega1_circular(dpotential,ddpotential,a)
        f2c = Omega2_circular(dpotential,a)

        # use the epicyclic approximation to step forward a tiny bit in radius
        f1h = Omega1_circular(dpotential,ddpotential,a+da)
        f2h = Omega2_circular(dpotential,a+da)

        # take a fairly large step in eccentricity to make sure we reach a safe zone
        # @ATTENTION, TOLECC needs to be smaller than de for this to work
        f1r,f2r = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,de,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        df1da = (f1h-f1c)/da
        df2da = (f2h-f2c)/da

        df1de = (f1r-f1c)/de
        df2de = (f2r-f2c)/de

        return f1c,f2c,df1da,df2da,df1de,df2de
end




"""compute_frequencies_rpra(ψ,dψ/dr,d²ψ/dr²,r_peri,r_apo[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e)
"""
function compute_frequencies_rpra(potential::Function,dpotential::Function,ddpotential::Function,
                                  r_peri::Float64,r_apo::Float64,TOLECC::Float64=0.001,verbose::Int64=0,NINT=32)

        f1,f2 = compute_frequencies_henon(potential,dpotential,ddpotential,r_peri,r_apo,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        return f1,f2
end
