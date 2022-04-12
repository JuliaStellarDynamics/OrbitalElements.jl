import OrbitalElements
using Printf

a = 1.0
e = 0.0

rp,ra = OrbitalElements.rpra_from_ae(a,e)
@printf("rp=%f ra=%f\n", rp,ra)

etaAE   = OrbitalElements.isochrone_eta_ae(rp,ra)
#print(etaAE)
f1real,f2real = OrbitalElements.isochrone_Omega_1_2(rp,ra)

f1comp,f2comp = OrbitalElements.compute_frequencies_henon_ae(OrbitalElements.isochrone_psi,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr,a,e)

@printf("O1=%f O1guess=%f O2=%f O2guess=%f\n", f1real,f1comp,f2real,f2comp)



#ae_from_omega1omega2_brute
E,L = OrbitalElements.isochrone_EL_from_alphabeta(f1real,f2real/f1real)
E,L = OrbitalElements.isochrone_EL_from_alphabeta(f1comp,f2comp/f1comp)

rpcomp,racomp = OrbitalElements.isochrone_rpra_fromEL(E,L)

@printf("rp=%f rpguess=%f ra=%f raguess=%f\n", rp,rpcomp,ra,racomp)

acomp,ecomp = OrbitalElements.ae_from_rpra(rpcomp,racomp)

acomp,ecomp = OrbitalElements.isochrone_ae_from_omega1omega2(f1comp,f2comp)
@printf("a=%f aguess=%f e=%f eguess=%f\n", a,acomp,e,ecomp)

aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1comp,f2comp,OrbitalElements.isochrone_psi,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr)

@printf("aguess=%f eguess=%f\n", aguess,eguess)



"""
aguess = OrbitalElements.Omega1circ_to_radius(f1comp,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr)
eguess = 0.5

for i=1:1000
    f1,f2,df1da,df2da,df1de,df2de = OrbitalElements.compute_frequencies_ae_derivs(OrbitalElements.isochrone_psi,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr,aguess,eguess)
    jacobian = [df1da df1de ; df2da df2de]
    omega1,omega2 = f1comp,f2comp
    step = jacobian \ (-([f1;f2] - [omega1 ; omega2]))
    aguess,eguess = aguess + step[1],eguess + step[2]
    @printf("aguess=%f eguess=%f\n", aguess,eguess)
end
"""
