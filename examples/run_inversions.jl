"""
check a bunch of quantities against the isochrone case, so we can be confident we are doing the numerical work correctly!
"""
import OrbitalElements
using Printf

# the (a,e) inversion breaks down between bc=8 and bc=10??
# how to redefine easy potentials to pass to frequency calculators
bc, M, G = 1.,1. ,1.
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

jrreal = OrbitalElements.isochrone_jr_rpra(rp,ra,bc,M,G)


f1comp,f2comp,a1comp = OrbitalElements.compute_frequencies_ae(potential,dpotential,ddpotential,a,e,true)


u = -1.
gval = OrbitalElements.Theta(potential,dpotential,ddpotential,u,rp,ra)

Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)
n1,n2 = 2,-3
rmax = 1000.
alpha,beta = f1real/Omega0,f2real/f1real
u,v = OrbitalElements.uv_from_alphabeta(alpha,beta,n1,n2,dpotential,ddpotential,rmax,Omega0)
alphaguess,betaguess = OrbitalElements.alphabeta_from_uv(u,v,n1,n2,dpotential,ddpotential,rmax,Omega0)

@printf("alpha=%f alphaguess=%f beta=%f betaguess=%f\n", alpha,alphaguess,beta,betaguess)

JacELab = OrbitalElements.isochrone_JacEL_to_alphabeta(alpha,beta)

# now try to get the empirical version of the jacobian

# get the estimated frequencies
f1rev,f2rev = alphaguess*Omega0,betaguess*alphaguess*Omega0
# get (a,e)
aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1rev,f2rev,potential,dpotential,ddpotential,0.000001,100)

# get the local derivs
f1,f2,df1da,df2da,df1de,df2de = OrbitalElements.compute_frequencies_ae_derivs(potential,dpotential,ddpotential,aguess,eguess)

# make the jacobian
#jaco1o2el = df1


# these should be rethought as change of variables. The Jacobian should be invertible, though?
# compute some nearby frequencies
f1c,f2c = OrbitalElements.compute_frequencies_henon_ae(potential,dpotential,ddpotential,aguess,eguess)
f1h,f2h = OrbitalElements.compute_frequencies_henon_ae(potential,dpotential,ddpotential,aguess+0.0001,eguess)
f1r,f2r = OrbitalElements.compute_frequencies_henon_ae(potential,dpotential,ddpotential,aguess,eguess+0.0001)
# and the corresponding E,L
rpc,rac = OrbitalElements.rpra_from_ae(aguess,eguess)
rph,rah = OrbitalElements.rpra_from_ae(aguess+0.0001,eguess)
rpr,rar = OrbitalElements.rpra_from_ae(aguess,eguess+0.0001)
Ec = OrbitalElements.E_from_rpra_pot(potential,dpotential,ddpotential,rpc,rac)
Lc = OrbitalElements.L_from_rpra_pot(potential,dpotential,ddpotential,rpc,rac)
Eh = OrbitalElements.E_from_rpra_pot(potential,dpotential,ddpotential,rph,rah)
Lr = OrbitalElements.L_from_rpra_pot(potential,dpotential,ddpotential,rpr,rar)
Lh = OrbitalElements.L_from_rpra_pot(potential,dpotential,ddpotential,rph,rah)

dEdf1 = (Eh-Ec)/(f1h-f1c)
dEdf2 = (Eh-Ec)/(f2h-f2c)
dLdf1 = (Lr-Lc)/(f1r-f1c)
dLdf2 = (Lr-Lc)/(f2r-f2c)
dLdf1 = (Lh-Lc)/(f1h-f1c)
dLdf2 = (Lh-Lc)/(f2h-f2c)

# this might be dangerous.
Jacguess = Omega0*f1c*(dEdf1*dLdf2 - dEdf2*dLdf1)


Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)
OrbitalElements.isochrone_EL_from_alphabeta(f1real/Omega0,f2real/f1real,bc,M,G)
rpcomp,racomp = OrbitalElements.isochrone_rpra_fromEL(E,L,bc,M,G)

@printf("rp=%f rpguess=%f ra=%f raguess=%f\n", rp,rpcomp,ra,racomp)

acomp,ecomp = OrbitalElements.ae_from_rpra(rpcomp,racomp)

acomp,ecomp = OrbitalElements.isochrone_ae_from_omega1omega2(f1comp,f2comp,bc,M,G)
@printf("a=%f aguess=%f e=%f eguess=%f\n", a,acomp,e,ecomp)

aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1comp,f2comp,potential,dpotential,ddpotential,0.000001,100)

@printf("aguess=%f eguess=%f\n", aguess,eguess)
