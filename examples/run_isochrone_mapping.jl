

"""
In contrast to the isochrone case, we do not have ground truth for Isochrone. So all comparisons will be made against high-resolution quantities.
"""


import OrbitalElements

# define easy potentials to pass to frequency calculators
bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)

println("Central pot=$(ψ(0.0))")
println("Central dpot=$(dψ(0.0))")
println("Central d2pot=$(d2ψ(0.0))")


rmin,rmax = 0.,10000.0
βc(αc::Float64)::Float64 = OrbitalElements.βcirc(αc,dψ,d2ψ,Ω₀,rmin=rmin,rmax=rmax)

αmin,αmax = OrbitalElements.αminmax(dψ,d2ψ,rmin,rmax,Ω₀=Ω₀)
println("(αmin,αmax)=($αmin,$αmax)")

n1,n2 = 10,-1
ωmin,ωmax = OrbitalElements.Findωminωmax(n1,n2,dψ,d2ψ,Ω₀=Ω₀,rmin=rmin,rmax=rmax)
println("(ωmin,ωmax)=($ωmin,$ωmax)")

uval = 0.5
vmin,vmax = OrbitalElements.FindVminVmax(uval,n1,n2,dψ,d2ψ,ωmin,ωmax,αmin,αmax,βc,Ω₀=Ω₀,rmin=rmin,rmax=rmax)
println("(vmin,vmax)=($vmin,$vmax)")

vval = 0.5*(vmax-vmin) + vmin


α,β = OrbitalElements.αβFromUV(uval,vval,n1,n2,ωmin,ωmax)
println("(u,v)=($uval,$vval)")
println("(α,β)=($α,$β)")

Ω1,Ω2 = α*Ω₀,α*β*Ω₀
# (Ω1,Ω2) -> (a,e)
a,e = OrbitalElements.AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ)
println("(a,e)=($a,$e)")


# select an (a,e) value for the orbit
a,e = 0.1, 0.1
