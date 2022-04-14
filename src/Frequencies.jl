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


"""compute_frequencies_ae(potential,dpotential,ddpotential,a,ecc[,TOLECC,verbose])

wrapper to select which type of frequency computation to perform, from (a,e)

"""
function compute_frequencies_ae(potential::Function,dpotential::Function,ddpotential::Function,
                                a::Float64,ecc::Float64,TOLECC::Float64=0.001,verbose::Int64=0)


        f1,f2 = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,ecc,TOLECC,verbose)

        return f1,f2
end


"""compute_ae_from_frequencies(potential,dpotential,ddpotential,a,ecc[,eps,maxiter])

wrapper to select which type of inversion to compute for (Omega1,Omega2)->(a,e)

"""
function compute_ae_from_frequencies(potential::Function,dpotential::Function,ddpotential::Function,
                                     omega1::Float64,omega2::Float64,eps::Float64=1*10^(-6),maxiter::Int64=10000)


        a,e = ae_from_omega1omega2_brute(omega1,omega2,potential,dpotential,ddpotential,eps,maxiter)

        return a,e
end


"""compute_frequencies_ae_derivs(potential,dpotential,ddpotential,a,ecc[,TOLECC,verbose])

wrapper to select which type of frequency computation to perform, from (a,e), but
DERIVATIVES

"""
function compute_frequencies_ae_derivs(potential::Function,
                                       dpotential::Function,
                                       ddpotential::Function,
                                       a::Float64,
                                       ecc::Float64,
                                       da::Float64=0.0001,
                                       de::Float64=0.0001,
                                       TOLECC::Float64=0.001,
                                       verbose::Int64=0)

        # grid is structured like
        # (f1h,f2h) [+da]
        #    ^
        # (f1c,f2c)-> (f1r,f2r) [+de]

        # @IMPROVE watch out for close to TOLECC, will fail across boundary
        f1c,f2c = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,ecc,TOLECC,verbose)

        f1h,f2h = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a+da,ecc,TOLECC,verbose)

        f1r,f2r = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,ecc+de,TOLECC,verbose)

        df1da = (f1h-f1c)/da
        df2da = (f2h-f2c)/da

        df1de = (f1r-f1c)/de
        df2de = (f2h-f2c)/de

        return f1c,f2c,df1da,df2da,df1de,df2de
end


"""compute_frequencies_rpra(potential,dpotential,ddpotential,r_peri,r_apo[,TOLECC,verbose])

wrapper to select which type of frequency computation to perform, from (a,e)

"""
function compute_frequencies_rpra(potential::Function,dpotential::Function,ddpotential::Function,
                                  r_peri::Float64,r_apo::Float64,TOLECC::Float64=0.001,verbose::Int64=0)

        f1,f2 = compute_frequencies_henon(potential,dpotential,ddpotential,r_peri,r_apo,TOLECC,verbose)

        return f1,f2
end
