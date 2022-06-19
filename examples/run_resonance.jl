"""
test resonance mapping capabilities
"""


import OrbitalElements
using Printf

# define easy potentials to pass to frequency calculators
const bc, M, G = 1.,1. ,1.
ψ       = r->OrbitalElements.isochrone_psi(r,bc,M,G)
dψdr    = r->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
d²ψdr²  = r->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Ω₀      = OrbitalElements.isochrone_Omega0(bc,M,G)

# this generates a function that computes beta (=Omega_2/Omega_1) as a function of Omega1
beta_c = OrbitalElements.make_betac(dψdr,d²ψdr²,2000,Ω₀)

# put in some dummy values for testing: picking a resonance
n1 = 2
n2 = 0
w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dψdr,d²ψdr²,1000.,Ω₀)
vbound = OrbitalElements.find_vbound(n1,n2,dψdr,d²ψdr²,1000.,Ω₀)

# for a given u value, find the integration boundaries
uval = -0.5
vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,beta_c)

beta_d(x::Float64) = x^(3/2)
vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,beta_d)
