"""
test some basic inversions to make sure all definitions are equivalent
"""


import OrbitalElements
using Printf

# define easy potentials to pass to frequency calculators
bc, M, G = 1.,1. ,1.
ψ       = r->OrbitalElements.isochrone_psi(r,bc,M,G)
dψdr    = r->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
d²ψdr²  = r->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Ω₀      = OrbitalElements.isochrone_Omega0(bc,M,G)


# select an (a,e) value for the orbit
a,e = 10., 0.4

# compute rperi and rapo
rp,ra = OrbitalElements.rpra_from_ae(a,e); @printf("rp=%f ra=%f\n", rp,ra)

Ω₁c,Ω₂c = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e)

# test the ability to recover (a,e) from the calculated frequencies
a1,e1 = OrbitalElements.compute_ae_from_frequencies(ψ,dψdr,d²ψdr²,Ω₁c,Ω₂c,1*10^(-12),1)

sma,ecc = OrbitalElements.compute_ae_from_frequencies(ψ,dψdr,d²ψdr²,Ω₁c,Ω₂c)

acirc = OrbitalElements.Omega1circ_to_radius(Ω₁c,dψdr,d²ψdr²)
# find that the problem is large a orbits, with low precision


E,L   = OrbitalElements.EL_from_rpra_pot(ψ,dψdr,d²ψdr²,rp,ra)
aguess,eguess = OrbitalElements.ae_from_EL_brute(E,L,ψ,dψdr,d²ψdr²,1*10^(-10),1000,0.001,0)

Ω₁c,Ω₂c = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,aguess,eguess)
aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(Ω₁c,Ω₂c,ψ,dψdr,d²ψdr²,0.000001,100)

# conversions from α, β to Omega1,Omega2
α,β = Ω₁c/Ω₀,Ω₂c/Ω₁c
Ω₁,Ω₂ = α*Ω₀,α*β*Ω₀


f1real,f2real = OrbitalElements.isochrone_Omega_1_2(rp,ra,bc,M,G)
Ω₁c,Ω₂c = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e)
@printf("O1=%f O1guess=%f O2=%f O2guess=%f\n", f1real,Ω₁c,f2real,Ω₂c)

# also compute actions?
jrreal = OrbitalElements.isochrone_jr_rpra(rp,ra,bc,M,G)
Ω₁c,Ω₂c,a1comp = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e,true)
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

# now try to get the empirical version of the jacobian

# get the estimated frequencies
f1rev,f2rev = αguess*Ω₀,βguess*αguess*Ω₀
# get (a,e)
aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1rev,f2rev,ψ,dψdr,d²ψdr²,0.000001,100)

# get estimates for local frequency derivatives the local derivs
f1,f2,df1da,df2da,df1de,df2de = OrbitalElements.compute_frequencies_ae_derivs(ψ,dψdr,d²ψdr²,aguess,eguess)
