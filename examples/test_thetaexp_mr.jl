using Plots
using LaTeXStrings

include("../src/Potentials/mestelzang.jl")

R0, V0 = 20., 1.
eps0 = 1e-1
ψ(r::Float64)::Float64   = mestel_psi(r,R0,V0,eps0)
dψ(r::Float64)::Float64  = mestel_dpsi_dr(r,R0,V0,eps0)
d2ψ(r::Float64)::Float64 = mestel_ddpsi_ddr(r,R0,V0,eps0)
d3ψ(r::Float64)::Float64 = mestel_d3psi_dr3(r,R0,V0,eps0)
d4ψ(r::Float64)::Float64 = mestel_d4psi_dr4(r,R0,V0,eps0)
# d4ψ(r::Float64)::Float64 = 0.0*r


include("../src/Utils/ComputeEL.jl")
include("../src/Henon/Frequencies.jl")
include("../src/Henon/Ufunc.jl")



# Array of results 
TOLECC=0.01
EDGE=0.02
a, e = 1.3, 0.45
umin, umax, δu = -1.0, -1.0+3*EDGE, 1e-5
tabu = collect(umin:δu:umax)
nu = length(tabu)

tabΘ = zeros(Float64,nu)

for i = 1:nu
    tabΘ[i] = ThetaAE(ψ,dψ,d2ψ,d3ψ,tabu[i],a,e;EDGE=EDGE,TOLECC=TOLECC)
end


pl=plot(tabu,tabΘ,title = "a = "*string(a)*";   e = "*string(e), lw = 2)
xlabel!(pl, L"$u$")
ylabel!(pl, L"$\Theta (u)$")
savefig(pl,"figures/Theta_a_"*string(a)*"_e_"*string(e)*".pdf")
