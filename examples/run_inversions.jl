import OrbitalElements
using Printf

# how to redefine easy potentials to pass to frequency calculators
bc, M, G = 8.,0.5 ,12.
potential   = r->OrbitalElements.isochrone_psi(r,bc,M,G)
dpotential  = r->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
ddpotential = r->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)

a,e = 5.0, 0.2

rp,ra = OrbitalElements.rpra_from_ae(a,e)
@printf("rp=%f ra=%f\n", rp,ra)

E  = OrbitalElements.E_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
Ei = OrbitalElements.isochrone_E_from_rpra(rp,ra,bc,M,G)

L  = OrbitalElements.L_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
Li = OrbitalElements.isochrone_L_from_rpra(rp,ra,bc,M,G)

f1real,f2real = OrbitalElements.isochrone_Omega_1_2(rp,ra,bc,M,G)
f1comp,f2comp = OrbitalElements.compute_frequencies_ae(potential,dpotential,ddpotential,a,e)
@printf("O1=%f O1guess=%f O2=%f O2guess=%f\n", f1real,f1comp,f2real,f2comp)

Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)
n1,n2 = 2,-3
rmax = 1000.
alpha,beta = f1real/Omega0,f2real/f1real
u,v = OrbitalElements.uv_from_alphabeta(alpha,beta,n1,n2,dpotential,ddpotential,rmax,Omega0)
alphaguess,betaguess = OrbitalElements.alphabeta_from_uv(u,v,n1,n2,dpotential,ddpotential,rmax,Omega0)

@printf("alpha=%f alphaguess=%f beta=%f betaguess=%f\n", alpha,alphaguess,beta,betaguess)




Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)
OrbitalElements.isochrone_EL_from_alphabeta(f1real/Omega0,f2real/f1real,bc,M,G)
rpcomp,racomp = OrbitalElements.isochrone_rpra_fromEL(E,L,bc,M,G)

@printf("rp=%f rpguess=%f ra=%f raguess=%f\n", rp,rpcomp,ra,racomp)

acomp,ecomp = OrbitalElements.ae_from_rpra(rpcomp,racomp)

acomp,ecomp = OrbitalElements.isochrone_ae_from_omega1omega2(f1comp,f2comp,bc,M,G)
@printf("a=%f aguess=%f e=%f eguess=%f\n", a,acomp,e,ecomp)

aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1comp,f2comp,potential,dpotential,ddpotential,0.000001,100)

@printf("aguess=%f eguess=%f\n", aguess,eguess)
