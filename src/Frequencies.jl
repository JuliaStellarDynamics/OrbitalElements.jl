#=

wrapper to select the methodology for computing the frequencies

only selecting Henon anomaly mapping for now,
but this is where one could select for different anomalies

=#

# bring in the frequency mapping
include("Henon/Frequencies.jl")
# bring in the anomaly mapping (i.e. f(u))
include("Henon/Ufunc.jl")

function compute_frequencies_ae(potential::Function,dpotential::Function,ddpotential::Function,
                                a::Float64,ecc::Float64,TOLECC::Float64=0.01,verbose::Int64=0)
        #=

        wrapper to select which type of frequency computation to perform, from (a,e)

        =#

        f1,f2 = compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,ecc,TOLECC,verbose)

        return f1,f2
end

function compute_frequencies_rpra(potential::Function,dpotential::Function,ddpotential::Function,
                                  r_peri::Float64,r_apo::Float64,TOLECC::Float64=0.01,verbose::Int64=0)
        #=

        wrapper to select which type of frequency computation to perform, from (rp,ra)

        =#

        f1,f2 = compute_frequencies_henon(potential,dpotential,ddpotential,r_peri,r_apo,TOLECC,verbose)

        return f1,f2
end
