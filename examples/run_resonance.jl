"""
test resonance mapping capabilities

julia run_resonance.jl
# 1.119036 seconds (2.02 M allocations: 113.069 MiB, 7.55% gc time, 99.99% compilation time)

julia --compile=min run_resonance.jl
# 0.008095 seconds (9.06 k allocations: 613.031 KiB)

julia --optimize=0 run_resonance.jl
# 1.006424 seconds (2.02 M allocations: 113.069 MiB, 12.97% gc time, 99.99% compilation time)

julia --optimize=0 --inline=no run_resonance.jl
# 1.972923 seconds (1.86 M allocations: 105.962 MiB, 20.93% gc time, 99.98% compilation time)

julia --compile=no run_resonance.jl
# 0.022154 seconds (9.06 k allocations: 613.031 KiB) # also a TON of complaints

julia --compile=min --inline=no run_resonance.jl

NOTES
-without --compile=min, there are _many_ allocations relating to finding vmin,vmax empirically (conversely, the analytic version is actually hurt by this!)

"""

import AstroBasis
import OrbitalElements
using Printf

# define easy potentials to pass to frequency calculators
const bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64       = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64    = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64  = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64  = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64  = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Ω₀      =    OrbitalElements.Omega0Isochrone(bc,M,G)

# this generates a function that computes beta (=Omega_2/Omega_1) as a function of Omega1
beta_c = OrbitalElements.make_betac(dψ,d2ψ,2000,Ω₀)
βc(alpha::Float64)::Float64 = OrbitalElements.beta_circ(alpha,dψ,d2ψ,Ω₀)

# put in some dummy values for testing: picking a resonance
n1 = -1
n2 = 2
@time w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dψ,d2ψ,1000.,Ω₀,Ziter=36)
@printf("wmin=%f wmax=%f\n", w_min,w_max)

@time vbound = OrbitalElements.find_vbound(n1,n2,dψ,d2ψ,1000.,Ω₀)
@printf("vbound=%f\n", vbound)

# for a given u value, find the integration boundaries
uval = -0.5

RH1,RH2 = OrbitalElements.RootOfHOmega(uval,w_min,w_max,n1,n2,vbound,beta_c)
@printf("RH1=%f RH2=%f\n", RH1,RH2)

C3 = OrbitalElements.ConstraintThree(uval,w_min,w_max,n1,n2)
@printf("C3=%f\n", C3)


@time vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,beta_c)
@printf("vmin=%f vmax=%f\n", vmin,vmax)
println(typeof(beta_c))

midv = (vmax+vmin)/2
vval = midv

#alpha,beta = OrbitalElements.alphabeta_from_uv(uval,midv,n1,n2,dψ,d2ψ)
alpha,beta = OrbitalElements.alphabeta_from_uv(uval,midv,n1,n2,w_min,w_max)

println("alpha=$alpha, beta=$beta")

omega1,omega2 = alpha*Ω₀,alpha*beta*Ω₀
a1,e1 = OrbitalElements.compute_ae_from_frequencies(ψ,dψ,d2ψ,omega1,omega2,1*10^(-12),1)
maxestep = 0.005
sma,ecc = OrbitalElements.compute_ae_from_frequencies(ψ,dψ,d2ψ,omega1,omega2,1*10^(-12),1000,0.001,0.0001,max(0.0001,0.001a1),min(max(0.0001,0.1a1*e1),maxestep),0)

rp,ra = OrbitalElements.rpra_from_ae(sma,ecc)
println("rp=$rp, ra=$ra")

gval = OrbitalElements.ThetaRpRa(ψ,dψ,d2ψ,uval,rp,ra,EDGE=0.02)
println("gval=$gval")


Lval = OrbitalElements.LFromRpRa(ψ,dψ,d2ψ,rp,ra)
Sigma, Delta = (ra+rp)*0.5, (ra-rp)*0.5

# Current location of the radius, r=r(u): isn't this exactly rp?
rval = Sigma + Delta*OrbitalElements.henon_f(uval)

# the velocity for integration
dt1du, dt2du = omega1*gval, (omega2 - Lval/(rval^(2)))*gval
println("dt1du=$dt1du, dt2du=$dt2du")

const rb = 10.
lmax,nmax = 2,10
lharmonic=2
basis = AstroBasis.CB73Basis_create(lmax=lmax, nmax=nmax,G=G,rb=rb)
AstroBasis.fill_prefactors!(basis)
AstroBasis.tabUl!(basis,lharmonic,rval)
println(basis.tabUl)

alpha,beta = OrbitalElements.alphabeta_from_uv(uval,vval,n1,n2,dψ,d2ψ)
println("alpha=$alpha, beta=$beta")

w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dψ,d2ψ,1000.,Ω₀)
Jacalphabeta = OrbitalElements.JacalphabetaToUV(n1,n2,w_min,w_max,vval) #(alpha,beta) -> (u,v)
JacEL        = OrbitalElements.IsochroneJacELtoAlphaBeta(alpha,beta)          #(E,L) -> (alpha,beta)
JacJ         = (1/omega1)                                #(J) -> (E,L)
dimensionl   = (1/Ω₀)                                # remove dimensionality
println("Jab=$Jacalphabeta JacEL=$JacEL, JacJ=$JacJ, dimensionl=$dimensionl")

function empirical_beta_c(x::Float64)::Float64
    return beta_c(x)
end

println(typeof(empirical_beta_c))

"""for the isochrone"""
function analytic_beta_c(x::Float64)::Float64
    return 1/(1 + x^(2/3))
end

@time vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,analytic_beta_c)
@printf("vmin=%f vmax=%f\n", vmin,vmax)
println(typeof(analytic_beta_c))

x = 0.5
@printf("Compare betac: empirical=%f analytical=%f\n", empirical_beta_c(x),analytic_beta_c(x))
