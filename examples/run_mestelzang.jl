"""
check a bunch of quantities against the isochrone case, so we can be confident we are doing the numerical work correctly!
"""


import OrbitalElements
using Printf

# define easy potentials to pass to frequency calculators
V0Mes = 1.    # Characteristic circular velocity``
R0Mes = 20. # Characteristic radii

epsMes = 0.01 # Potential taper in the center

# RinMes, RoutMes, R0Mes, RmaxMes = 1., 11.5, 20., 20. # Characteristic radii
# xiMes = 0.5   # Active fraction
# qMes  = 11.44 # Power index
# nuMes, muMes = 4, 5 # Inner and outer tapering exponents


ψ(r::Float64)       = OrbitalElements.mestel_psi(r,R0Mes,V0Mes,epsMes)
dψdr(r::Float64)    = OrbitalElements.mestel_dpsi_dr(r,R0Mes,V0Mes,epsMes)
d²ψdr²(r::Float64)  = OrbitalElements.mestel_ddpsi_ddr(r,R0Mes,V0Mes,epsMes)
Ω₀      = OrbitalElements.mestel_Omega0(R0Mes,V0Mes,epsMes)

println(" Ω₀ : ",Ω₀)


# select an (a,e) value for the orbit
a,e =0.008, 0.4

# compute rperi and rapo
println("Compute rp,ra...")
@time rp,ra = OrbitalElements.rpra_from_ae(a,e); @printf("rp=%f ra=%f\n", rp,ra)

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
@time aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(Ω₁c,Ω₂c,ψ,dψdr,d²ψdr²,1*10^(-12),1000)
println("...",aguess,",",eguess)


# test the ability to recover (a,e) from the calculated frequencies
println("Compute a1,e1 (inversion)...")
#@time a1,e1 = OrbitalElements.compute_ae_from_frequencies(ψ,dψdr,d²ψdr²,Ω₁c,Ω₂c,1*10^(-12),100)
@time a1,e1 = OrbitalElements.compute_ae_from_frequencies(ψ,dψdr,d²ψdr²,Ω₁c,Ω₂c,1*10^(-12),1000)
println("...",a1,",",e1)


println("Compute sma,ecc...")
@time sma,ecc = OrbitalElements.compute_ae_from_frequencies(ψ,dψdr,d²ψdr²,Ω₁c,Ω₂c)

println("Compute acirc...")
@time acirc = OrbitalElements.Omega1circ_to_radius(Ω₁c,dψdr,d²ψdr²)
# find that the problem is large a orbits, with low precision

#E,L   = OrbitalElements.EL_from_rpra_pot(ψ,dψdr,d²ψdr²,rp,ra)
println("INVERT E,L...")
@time aguess,eguess = OrbitalElements.ae_from_EL_brute(E,L,ψ,dψdr,d²ψdr²,1*10^(-6),1000,0.001,0)

@printf("aguess=%f eguess=%f\n", aguess,eguess)


# conversions from α, β to Omega1,Omega2
α,β = Ω₁c/Ω₀,Ω₂c/Ω₁c
Ω₁,Ω₂ = α*Ω₀,α*β*Ω₀

Ω₁c,Ω₂c = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e)
@printf("O1guess=%f O2guess=%f\n", Ω₁c,Ω₂c)

# also compute actions?
Ω₁c,Ω₂c,a1comp = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e,action=true)
@printf("Jrguess=%f\n", a1comp)





u = -1.
v = 0.5
gval = OrbitalElements.Theta(ψ,dψdr,d²ψdr²,u,rp,ra)

n1,n2 = 2,-3
rmax = 1000.
αguess,βguess = OrbitalElements.alphabeta_from_uv(u,v,n1,n2,dψdr,d²ψdr²,rmax,Ω₀)
uguess,vguess = OrbitalElements.uv_from_alphabeta(αguess,βguess,n1,n2,dψdr,d²ψdr²,rmax,Ω₀)

@printf("αguess=%f βguess=%f\n", αguess,βguess)
@printf("u=%f uguess=%f v=%f vguess=%f\n", u,uguess,v,vguess)

JacELab = OrbitalElements.isochrone_JacEL_to_alphabeta(α,β)

# now try to get the empirical version of the jacobian

# get the estimated frequencies
f1rev,f2rev = αguess*Ω₀,βguess*αguess*Ω₀
# get (a,e)
aguess,eguess = OrbitalElements.ae_from_omega1omega2_brute(f1rev,f2rev,ψ,dψdr,d²ψdr²,0.000001,100)

@printf("aguess=%f eguess=%f\n", aguess,eguess)

# get estimates for local frequency derivatives the local derivs
f1,f2,df1da,df2da,df1de,df2de = OrbitalElements.compute_frequencies_ae_derivs(ψ,dψdr,d²ψdr²,aguess,eguess)



