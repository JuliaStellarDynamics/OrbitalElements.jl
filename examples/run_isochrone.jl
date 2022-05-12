"""
check a bunch of quantities against the isochrone case, so we can be confident we are doing the numerical work correctly!
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
a,e = 1., 0.4

# compute rperi and rapo
rp,ra = OrbitalElements.rpra_from_ae(a,e); @printf("rp=%f ra=%f\n", rp,ra)

# test frequency computation
Ω₁r,Ω₂r = OrbitalElements.isochrone_Omega_1_2(rp,ra,bc,M,G)
Ω₁c,Ω₂c = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e)
@printf("O1=%f O1guess=%f O2=%f O2guess=%f\n", Ω₁r,Ω₁c,Ω₂r,Ω₂c)

# compute (alpha,beta)
alpha,beta = Ω₁r/Ω₀,Ω₂r/Ω₁r

# invert to get E,L: this can still possibly fail for bad values of sqrts.
E,L = OrbitalElements.isochrone_EL_from_alphabeta(alpha,beta,bc,M,G)
rpcomp,racomp = OrbitalElements.isochrone_rpra_fromEL(E,L,bc,M,G)

# compare to original values
@printf("rp=%f rpguess=%f ra=%f raguess=%f\n", rp,rpcomp,ra,racomp)
acomp,ecomp = OrbitalElements.ae_from_rpra(rpcomp,racomp)
@printf("empircal: a=%f aguess=%f e=%f eguess=%f\n", a,acomp,e,ecomp)

acomp,ecomp = OrbitalElements.isochrone_ae_from_omega1omega2(Ω₁c,Ω₂c,bc,M,G)
@printf("analytical: a=%f aguess=%f e=%f eguess=%f\n", a,acomp,e,ecomp)

# more options for checking...
#E,L   = OrbitalElements.EL_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
#Ei,Li = OrbitalElements.isochrone_EL_from_rpra(rp,ra,bc,M,G)
