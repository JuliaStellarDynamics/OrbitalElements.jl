"""
test some basic inversions to make sure all definitions are equivalent

julia --compile=min run_inversions.jl

"""

import OrbitalElements
using Printf

# define easy potentials to pass to frequency calculators
#bc, M, G = 1.,1. ,1.
const bc, M, G = 1.,1. ,1.  # these can be constant for optimization flags

Ω₀                  = OrbitalElements.isochrone_Omega0(bc,M,G)

#ψ(r::Float64)       = r->OrbitalElements.isochrone_psi(r,bc,M,G)
#dψdr(r::Float64)    = r->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
#d²ψdr²(r::Float64)  = r->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)

#ψ       = r->OrbitalElements.isochrone_psi(r,bc,M,G)
#dψdr    = r->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
#d²ψdr²  = r->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)

#ψ       = r::Float64->OrbitalElements.isochrone_psi(r,bc,M,G)
#dψdr    = r::Float64->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
#d²ψdr²  = r::Float64->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)

ψ(r::Float64)::Float64       = OrbitalElements.isochrone_psi(r,bc,M,G)
dψdr(r::Float64)::Float64    = OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
d²ψdr²(r::Float64)::Float64  = OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)

# select an (a,e) value for the orbit
a,e = .1, 0.4

# compute rperi and rapo
println("Compute rp,ra...")
@time rp,ra = OrbitalElements.rpra_from_ae(a,e); @printf("rp=%f ra=%f\n", rp,ra)


# compute circular equivalent
println("Compute O1 circular...")
O1c   = OrbitalElements.Omega1_circular(dψdr,d²ψdr²,ra)
println(O1c)

# compute (E,L)
println("Compute E,L...")
E,L   = OrbitalElements.EL_from_rpra_pot(ψ,dψdr,d²ψdr²,rp,ra)

# compute frequencies
println("Compute Ω₁c,Ω₂c...")
@time Ω₁c,Ω₂c = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e)
println("...",Ω₁c,",",Ω₂c)
println("Compute Ω₁c,Ω₂c (Hénon specific)...")
@time Ω₁c,Ω₂c = OrbitalElements.henon_anomaly_frequencies(ψ,ra,rp,E,L,NINT=64)
println("...",Ω₁c,",",Ω₂c)
@time Ω₁c,Ω₂c = OrbitalElements.henon_anomaly_frequencies(ψ,ra,rp,E,L,NINT=32)
println("...",Ω₁c,",",Ω₂c)
@time Ω₁c,Ω₂c = OrbitalElements.henon_anomaly_frequencies(ψ,ra,rp,E,L,NINT=16)
println("...",Ω₁c,",",Ω₂c)

# do the frequency inversion.
Ω₁c,Ω₂c = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e)
println("INVERT Ω₁c,Ω₂c...")
@time aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(Ω₁c,Ω₂c,ψ,dψdr,d²ψdr²,1*10^(-12),100)
println("...",aguess,",",eguess)


# test the ability to recover (a,e) from the calculated frequencies
println("Compute a1,e1 (inversion)...")
#@time a1,e1 = OrbitalElements.compute_ae_from_frequencies(ψ,dψdr,d²ψdr²,Ω₁c,Ω₂c,1*10^(-12),100)
@time a1,e1 = OrbitalElements.compute_ae_from_frequencies(ψ,dψdr,d²ψdr²,Ω₁c,Ω₂c,1*10^(-12),100)
println("...",a1,",",e1)


Ω₁c,Ω₂c = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e)
println("INVERT Ω₁c,Ω₂c...part 2")
@time aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(Ω₁c,Ω₂c,ψ,dψdr,d²ψdr²,1*10^(-12),100)
println("...",aguess,",",eguess)

println("Compute sma,ecc...")
@time sma,ecc = OrbitalElements.compute_ae_from_frequencies(ψ,dψdr,d²ψdr²,Ω₁c,Ω₂c)
println("...",sma,",",ecc)

println("Compute acirc...")
@time acirc = OrbitalElements.Omega1circ_to_radius(Ω₁c,dψdr,d²ψdr²)
# find that the problem is large a orbits, with low precision

#E,L   = OrbitalElements.EL_from_rpra_pot(ψ,dψdr,d²ψdr²,rp,ra)
println("INVERT E,L...")
@time aguess,eguess = OrbitalElements.ae_from_EL_brute(E,L,ψ,dψdr,d²ψdr²,1*10^(-10),1000,0.001,0)




# conversions from α, β to Omega1,Omega2
α,β = Ω₁c/Ω₀,Ω₂c/Ω₁c
Ω₁,Ω₂ = α*Ω₀,α*β*Ω₀


f1real,f2real = OrbitalElements.IsochroneOmega12FromRpRa(rp,ra,bc,M,G)
Ω₁c,Ω₂c = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e)
@printf("O1=%f O1guess=%f O2=%f O2guess=%f\n", f1real,Ω₁c,f2real,Ω₂c)

# also compute actions?
jrreal = OrbitalElements.isochrone_jr_rpra(rp,ra,bc,M,G)
Ω₁c,Ω₂c,a1comp = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e,action=true)
@printf("Jr=%f Jrguess=%f\n", jrreal,a1comp)





u = -1.
gval = OrbitalElements.Theta(ψ,dψdr,d²ψdr²,u,rp,ra)

Ω₀ = OrbitalElements.isochrone_Omega0(bc,M,G)
n1,n2 = 2,-3
rmax = 1000.
α,β = f1real/Ω₀,f2real/f1real
u,v = OrbitalElements.uv_from_alphabeta(α,β,n1,n2,dψdr,d²ψdr²,rmax,Ω₀)
αguess,βguess = OrbitalElements.alphabeta_from_uv(u,v,n1,n2,dψdr,d²ψdr²,rmax,Ω₀)

@printf("α=%f αguess=%f β=%f βguess=%f\n", α,αguess,β,βguess)

JacELab = OrbitalElements.isochrone_JacEL_to_alphabeta(α,β)
JacELab2 = OrbitalElements.JacELToAlphaBetaAE(a,e,ψ,dψdr,d²ψdr²,Ω₀)
println("Compare Jacobians:JacELab=$JacELab,JacELab2=$JacELab2")
# now try to get the empirical version of the jacobian
