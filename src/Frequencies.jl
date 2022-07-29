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

@ATTENTION can use the isochrone-specific if you are using an isochrone. Otherwise this is a bit costly.
"""
function JacEL_to_alphabeta(alpha::Float64,beta::Float64)
    isochrone_JacEL_to_alphabeta(alpha,beta)
end

"""

@ATTENTION this takes (a,e) as arguments.
@ATTENTION this combines several numerical derivatives; please take care!

@IMPROVE add massaging parameters for numerical derivatives
@IMPROVE fix boundary values when using limited development
@IMPROVE noisy at the boundaries
"""
function JacELToAlphaBetaAE(a::Float64,ecc::Float64,ψ::Function,dψdr::Function,d²ψdr²::Function,Ω₀::Float64)

    # to be fixed for limited development...
    if ecc>0.99
        ecc=0.99
    end

    if ecc<0.01
        ecc=0.01
    end

    # get all numerical derivatives
    f1c,f2c,df1da,df2da,df1de,df2de = compute_frequencies_ae_derivs(ψ,dψdr,d²ψdr²,a,ecc)
    Ec,Lc,dEda,dEde,dLda,dLde = dEdL_from_ae_pot(ψ,dψdr,d²ψdr²,a,ecc)

    # construct Jacobians
    J_EL_ae = abs(dEda*dLde - dEde*dLda)
    J_o1o2_ae = abs(df1da*df2de - df1de*df2da)

    # combine and return
    return f1c*Ω₀*J_EL_ae/J_o1o2_ae

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


"""compute_frequencies_ae_derivs(ψ,dψ/dr,d²ψ/dr²a,ecc[,TOLECC,verbose])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES

@IMPROVE: could add action derivatives here? more copacetic with analytic derivatives anyway
"""
function compute_frequencies_ae_derivs(potential::Function,
                                       dpotential::Function,
                                       ddpotential::Function,
                                       a::Float64,
                                       ecc::Float64,
                                       da::Float64=0.0001,
                                       de::Float64=0.0001,
                                       TOLECC::Float64=0.001,
                                       verbose::Int64=0,
                                       NINT::Int64=32)

        # grid is structured like
        # (f1h,f2h) [+da]
        #    ^
        # (f1c,f2c)-> (f1r,f2r) [+de]

        # @IMPROVE watch out for close to TOLECC, will fail across boundary
        f1c,f2c = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,ecc,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        f1h,f2h = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a+da,ecc,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

        f1r,f2r = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,ecc+de,action=false,TOLECC=TOLECC,verbose=verbose,NINT=NINT)

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
