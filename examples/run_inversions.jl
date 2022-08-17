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
#a,e = 2.6636008542719694,0.9824936963436619
#a,e = 100.,0.001
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

# if you want to try inverting a specific frequency set, put them here:
#Ω₁a,Ω₂a = 0.1757002803595819,0.08954145469595597

# try an inversion back to (a,e) using third derivative
@time aa,ea = OrbitalElements.AEFromOmega1Omega2Brute(Ω₁a,Ω₂a,ψ,dψ,d2ψ,d3ψ,NINT=32,EDGE=0.01,verbose=3,ITERMAX=10)
println("Recovered a=$aa e=$ea")

@time al,el = OrbitalElements.ae_from_omega1omega2_brute(Ω₁a,Ω₂a,ψ,dψ,d2ψ)
println("LazyR     a=$al e=$el")
