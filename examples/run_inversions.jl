import OrbitalElements
using Printf

bc, M, G = 2., 1. ,1.
potential   = r->OrbitalElements.isochrone_psi(r,bc,M,G)
dpotential  = r->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
ddpotential = r->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)

a,e = 2.0, 0.5

rp,ra = OrbitalElements.rpra_from_ae(a,e)
@printf("rp=%f ra=%f\n", rp,ra)


f1real,f2real = OrbitalElements.isochrone_Omega_1_2(rp,ra,bc,M,G)
f1comp,f2comp = OrbitalElements.compute_frequencies_henon_ae(potential,dpotential,ddpotential,a,e)
@printf("O1=%f O1guess=%f O2=%f O2guess=%f\n", f1real,f1comp,f2real,f2comp)

# M and G scaling are fine, bc is not.


E,L = OrbitalElements.isochrone_EL_from_alphabeta(f1real,f2real/f1real)
E,L = OrbitalElements.isochrone_EL_from_alphabeta(f1comp,f2comp/f1comp)

rpcomp,racomp = OrbitalElements.isochrone_rpra_fromEL(E,L)

@printf("rp=%f rpguess=%f ra=%f raguess=%f\n", rp,rpcomp,ra,racomp)

acomp,ecomp = OrbitalElements.ae_from_rpra(rpcomp,racomp)

acomp,ecomp = OrbitalElements.isochrone_ae_from_omega1omega2(f1comp,f2comp)
@printf("a=%f aguess=%f e=%f eguess=%f\n", a,acomp,e,ecomp)

aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1comp,f2comp,OrbitalElements.isochrone_psi,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr)

@printf("aguess=%f eguess=%f\n", aguess,eguess)
