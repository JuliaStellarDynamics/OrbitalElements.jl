"""
test resonance mapping capabilities

"""

import OrbitalElements

# define easy potentials to pass to frequency calculators
const bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64    = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64   = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64  = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64  = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64  = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Ω₀      =    OrbitalElements.Omega0Isochrone(bc,M,G)

function check_mapping(u::Float64,v::Float64,n1::Int64,n2::Int64)
    w_min,w_max = OrbitalElements.FindWminWmax(n1,n2,dψ,d2ψ,1000.,Ω₀,Ziter=36)
    println("Empirical wmin=$w_min,wmax=$w_max")

    w_min,w_max = OrbitalElements.FindWminWmaxIsochrone(n1,n2)
    println("True      wmin=$w_min,wmax=$w_max")

    alpha,beta = OrbitalElements.AlphaBetaFromUV(u,v,n1,n2,w_min,w_max)
    println("Alpha=$alpha, Beta=$beta")

    omega1,omega2 = alpha*Ω₀,alpha*beta*Ω₀
    sma,ecc = OrbitalElements.IsochroneAEFromOmega1Omega2(omega1,omega2,bc,M,G)
    println("True      a=$sma,e=$ecc")

    sma,ecc = OrbitalElements.AEFromOmega1Omega2Brute(omega1,omega2,ψ,dψ,d2ψ,d3ψ,NINT=32,EDGE=0.01,verbose=0)
    println("Empirical a=$sma,e=$ecc")
end

u = -0.9999
n1,n2 = -1,2
vmin,vmax = OrbitalElements.FindVminVmaxIsochrone(n1,n2,u)
println("Analytic  vmin=$vmin,vmax=$vmax")

ωmin,ωmax = OrbitalElements.FindWminWmax(n1,n2,dψ,d2ψ,10000.,Ω₀)
beta_c(alpha_c::Float64)::Float64 = OrbitalElements.beta_circ(alpha_c,dψ,d2ψ,Ω₀,rmax=1.0e6)
vbound = OrbitalElements.FindVbound(n1,n2,dψ,d2ψ,10000.,Ω₀)
#vmin,vmax = OrbitalElements.FindVminVmax(u,ωmin,ωmax,n1,n2,vbound,beta_c)
println("Empirical vmin=$vmin,vmax=$vmax")

check_mapping(u,(vmin+vmax)/2,n1,n2)

#=
vmin,vmax = OrbitalElements.FindVminVmaxIsochrone(3,n2,u)
println("vmin=$vmin,vmax=$vmax")
check_mapping(u,(vmin+vmax)/2,3,n2)

vmin,vmax = OrbitalElements.FindVminVmaxIsochrone(4,n2,u)
println("vmin=$vmin,vmax=$vmax")
check_mapping(u,(vmin+vmax)/2,4,n2)
=#
