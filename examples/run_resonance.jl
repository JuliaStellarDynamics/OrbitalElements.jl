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


import OrbitalElements
using Printf

# define easy potentials to pass to frequency calculators
const bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64       = OrbitalElements.isochrone_psi(r,bc,M,G)
dψdr(r::Float64)::Float64    = OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
d²ψdr²(r::Float64)::Float64  = OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Ω₀      =    OrbitalElements.isochrone_Omega0(bc,M,G)

# this generates a function that computes beta (=Omega_2/Omega_1) as a function of Omega1
beta_c = OrbitalElements.make_betac(dψdr,d²ψdr²,2000,Ω₀)

# put in some dummy values for testing: picking a resonance
n1 = 2
n2 = 0
@time w_min,w_max = OrbitalElements.find_wmin_wmax(n1,n2,dψdr,d²ψdr²,1000.,Ω₀)
@printf("wmin=%f wmax=%f\n", w_min,w_max)

@time vbound = OrbitalElements.find_vbound(n1,n2,dψdr,d²ψdr²,1000.,Ω₀)
@printf("vbound=%f\n", vbound)

# for a given u value, find the integration boundaries
uval = -0.5
@time vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,beta_c)
@printf("vmin=%f vmax=%f\n", vmin,vmax)
println(typeof(beta_c))

function empirical_beta_c(x::Float64)::Float64
    return beta_c(x)
end

println(typeof(empirical_beta_c))


function analytic_beta_c(x::Float64)::Float64
    return 1/(1 + x^(2/3))
end

@time vmin,vmax = OrbitalElements.find_vmin_vmax(uval,w_min,w_max,n1,n2,vbound,analytic_beta_c)
@printf("vmin=%f vmax=%f\n", vmin,vmax)
println(typeof(analytic_beta_c))

x = 0.5
@printf("Compare betac: empirical=%f analytical=%f\n", empirical_beta_c(x),analytic_beta_c(x))
