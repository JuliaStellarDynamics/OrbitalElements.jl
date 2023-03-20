
import OrbitalElements
using BenchmarkTools


######################################################################
#
# Potential to consider
#
######################################################################

#####
# Isochrone
#####
const bc, M, G = 1.,1. ,1.
const ψ(r::Float64)    = OrbitalElements.ψIsochrone(r,bc,M,G)
const dψ(r::Float64)   = OrbitalElements.dψIsochrone(r,bc,M,G)
const d2ψ(r::Float64)  = OrbitalElements.d2ψIsochrone(r,bc,M,G)
const d3ψ(r::Float64)  = OrbitalElements.d3ψIsochrone(r,bc,M,G)
const d4ψ(r::Float64)  = OrbitalElements.d4ψIsochrone(r,bc,M,G)
const Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)

const params = OrbitalElements.OrbitalParameters()


######################################################################
#
# Variables values
#
######################################################################

n1, n2 = -1, 2
a, e = 1.0, 0.05
u = 0.4
Ω1, Ω2 = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
E = OrbitalElements.EFromAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
J, L  = OrbitalElements.ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
rp, ra = OrbitalElements.RpRaFromAE(a,e)
r = OrbitalElements.ru(u,a,e)

println("UNDEFINED PARAMETER STRUCTURE BENCHMARKS:")
println("Frequencies integrand computation Benchmark")
@btime OrbitalElements.ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e)

println("Frequencies computation Benchmark")
@btime OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)

println("Frequencies inversion Benchmark")
@btime OrbitalElements.ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,d4ψ,Ω1,Ω2)

#####
# Defined parameter structure
#####
println("DEFINED PARAMETER STRUCTURE BENCHMARKS:")
println("Frequencies integrand computation Benchmark")
@btime OrbitalElements.ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,params)

println("Frequencies computation Benchmark")
@btime OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)

println("Frequencies inversion Benchmark")
@btime OrbitalElements.ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,d4ψ,Ω1,Ω2,params)
######################################################################
#
# Benchmarks (in REPL)
#
######################################################################

# Open REPL at the directory location:
#   - (base) me@mycomputer OrbitalElements/ % julia
# Activate the environment of the package:
#   - julia> ]                              # Open package manager (pkg)
#   - (v1.7) pkg> activate .                # Activate the local environment
#   - (OrbitalElements) pkg> (backspace)    # Close package manager
#   - julia>
# Include the test file:
#   - julia> include("test/test.jl")
#
#
# Perform the following possible benchmarks:

#####
# Frequency integrand
#####
# @benchmark OrbitalElements.ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,params)                                # (dans src/Henon/Ufunc.jl)

# Expected result:
# Memory estimate: 0 bytes, allocs estimate: 0

#####
# Frequency computation
#####
# @benchmark OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)             # (dans src/Henon/Frequencies.jl)
# Equivalent à (appelle)
# @benchmark OrbitalElements.αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)              # (dans src/Frequencies.jl)

# Expected result:
# Memory estimate: 0 bytes, allocs estimate: 0

#####
# Frequency and dericatives computation
#####
# @benchmark OrbitalElements.ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
# (dans src/Frequencies.jl)        # Dérivées à l'extérieur (the one used in pratice) 
# Pas vraiment équivalent à 
# @benchmark OrbitalElements.DFrequenciesHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
# (dans src/Henon/Frequencies.jl)   # Derivées sous l'intégrale            


# Expected result:
# ComputeFrequenciesAEWithDeriv -> Memory estimate: 0 bytes, allocs estimate: 0.
# DFrequenciesHenonΘAE          -> Memory estimate: 0 bytes, allocs estimate: 0.


#####
# Inversion (Ω1,Ω2) -> (a,e)
#####
# @benchmark OrbitalElements.AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ,d4ψ,params)                
# (dans src/Utils/NumericalInversion.jl)
# équivalent à (est appelé par)
# @benchmark OrbitalElements.ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,d4ψ,Ω1,Ω2,params)       
# (dans src/Frequencies.jl)

# Expected result:
# AEFromΩ1Ω2Brute           -> Memory estimate: 0 bytes, allocs estimate: 0.
# ComputeAEFromFrequencies  -> Memory estimate: 0 bytes, allocs estimate: 0.