"""
test some basic inversions to make sure all definitions are equivalent

julia --compile=min run_inversions.jl

"""

using BenchmarkTools
import OrbitalElements

# define easy potentials to pass to frequency calculators
#bc, M, G = 1.,1. ,1.
const bc, M, G = 1.,1. ,1.  # these can be constant for optimization flags

ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
Ω₀                       = OrbitalElements.Omega0Isochrone(bc,M,G)

# select an (a,e) value for the orbit
a,e = 0.005, 0.9
println("Input     a=$a e=$e")


# compute rperi and rapo
#println("Compute rp,ra...")
#rp,ra = OrbitalElements.rpra_from_ae(a,e)
#println("rp=$rp ra=$ra")

# compute exact frequencies
Ω₁e,Ω₂e = OrbitalElements.IsochroneOmega12FromAE(a,e,bc,M,G)
println("Exact       Ω₁=$Ω₁e, Ω₂=$Ω₂e")

# compute approximate frequencies
@time Ω₁a,Ω₂a,Jrr = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e,NINT=512,action=true)
println("Approximate Ω₁=$Ω₁a, Ω₂=$Ω₂a")

# try an inversion back to (a,e) using third derivative
@time aa,ea = OrbitalElements.AEFromOmega1Omega2Brute(Ω₁a,Ω₂a,ψ,dψ,d2ψ,d3ψ,NINT=32,EDGE=0.01,verbose=1)
println("Recovered a=$aa e=$ea")

@time al,el = OrbitalElements.ae_from_omega1omega2_brute(Ω₁a,Ω₂a,ψ,dψ,d2ψ)
println("LazyR     a=$al e=$el")
