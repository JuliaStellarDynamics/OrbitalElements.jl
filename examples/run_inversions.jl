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

a,e = 10., 0.4
rp,ra = OrbitalElements.rpra_from_ae(a,e); @printf("rp=%f ra=%f\n", rp,ra)

f1comp,f2comp = OrbitalElements.compute_frequencies_ae(potential,dpotential,ddpotential,a,e) # frequencies go over the limit...hopefully not a problem
a1,e1 = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,f1comp,f2comp,1*10^(-12),1)

maxestep = 0.005
aguess,eguess = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,f1comp,f2comp,1*10^(-12),1000,0.001,0.0001,max(0.0001,0.001a1),min(max(0.0001,0.1a1*e1),maxestep),1)


aguess,eguess = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,f1comp,f2comp,1*10^(-12),5000,0.001,0.0001,max(0.0001,0.001a1),min(max(0.0001,0.1a1*e1),maxestep),1)
aguess,eguess = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,f1comp,f2comp,1*10^(-12),5000,0.001,0.0001,max(0.0001,0.005a1),min(max(0.0001,0.1a1*e1),maxestep),1)
aguess,eguess = OrbitalElements.compute_ae_from_frequencies(potential,dpotential,ddpotential,f1comp,f2comp,1*10^(-12),5000,0.001,0.0001,max(0.0001,0.01a),min(max(0.0001,0.1a*e),maxestep),1)


aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1comp,f2comp,potential,dpotential,ddpotential,1*10^(-12),5000,0.001,0.0001,max(0.0001,0.01a),min(max(0.0001,0.1a*e),0.01),0)


# works normally
aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1comp,f2comp,potential,dpotential,ddpotential,1*10^(-12),5000,0.001,0.0001,max(0.0001,0.001a),max(0.0001,0.01a*e),0)



E,L   = OrbitalElements.EL_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
Ei,Li = OrbitalElements.isochrone_EL_from_rpra(rp,ra,bc,M,G)

aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1comp,f2comp,potential,dpotential,ddpotential,1*10^(-8),500,0.001,0.0001,max(0.0001,0.001a),0.01e,0)


da,de=0.0001,0.0001
da,de=0.001,0.001
a,e = 50., 0.95

OrbitalElements.compute_frequencies_ae_derivs(potential,dpotential,ddpotential,a,e,da,de,0.001,0)


aguess,eguess = OrbitalElements.ae_from_EL_brute(E,L,potential,dpotential,ddpotential,1*10^(-10),1000,0.001,0)

f1comp,f2comp = OrbitalElements.compute_frequencies_ae(potential,dpotential,ddpotential,aguess,eguess)
aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1comp,f2comp,potential,dpotential,ddpotential,0.000001,100)


# test out a numerical derivative
aguessh,eguessh = OrbitalElements.ae_from_EL_brute(E+0.0001,L,potential,dpotential,ddpotential,1*10^(-10),1000,0.001,0)
aguessr,eguessr = OrbitalElements.ae_from_EL_brute(E,L+0.0001,potential,dpotential,ddpotential,1*10^(-10),1000,0.001,0)
dadE = (aguessh-aguess)/0.0001
dedE = (eguessh-eguess)/0.0001
dadL = (aguessr-aguess)/0.0001
dedL = (eguessr-eguess)/0.0001

# get the local derivs
f1,f2,df1da,df2da,df1de,df2de = OrbitalElements.compute_frequencies_ae_derivs(potential,dpotential,ddpotential,aguess,eguess)

df1dE = dadE*df1da + dedE*df1de
df1dL = dadL*df1da + dedL*df1de
df2dE = dadE*df2da + dedE*df2de
df2dL = dadL*df2da + dedL*df2de

Jac_f1f2_EL = abs(df1dE*df2dL - df1dL*df2dE)


Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)
alpha,beta = f1/Omega0,f2/f1
JacELab = OrbitalElements.isochrone_JacEL_to_alphabeta(alpha,beta)

estJacELab = Omega0*f1*Jac_f1f2_EL

Eguess,Lguess,dEda,dEde,dLda,dLde = OrbitalElements.dEdL_from_rpra_pot(potential,dpotential,ddpotential,rp,ra,0.0001,0.0001,0.001)
df1dE = df1da/dEda + df1de/dEde
df1dL = df1da/dLda + df1de/dLde
df2dE = df2da/dEda + df2de/dEde
df2dL = df2da/dLda + df2de/dLde


Jac_f1f2_EL = abs(df1dE*df2dL - df1dL*df2dE)


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



Omega0 = OrbitalElements.isochrone_Omega0(bc,M,G)
OrbitalElements.isochrone_EL_from_alphabeta(f1real/Omega0,f2real/f1real,bc,M,G)
rpcomp,racomp = OrbitalElements.isochrone_rpra_fromEL(E,L,bc,M,G)

@printf("rp=%f rpguess=%f ra=%f raguess=%f\n", rp,rpcomp,ra,racomp)

acomp,ecomp = OrbitalElements.ae_from_rpra(rpcomp,racomp)

acomp,ecomp = OrbitalElements.isochrone_ae_from_omega1omega2(f1comp,f2comp,bc,M,G)
@printf("a=%f aguess=%f e=%f eguess=%f\n", a,acomp,e,ecomp)

aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1comp,f2comp,potential,dpotential,ddpotential,0.000001,100)

@printf("aguess=%f eguess=%f\n", aguess,eguess)
