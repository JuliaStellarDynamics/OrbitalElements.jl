version = "new"

using HDF5

include("../src/Potentials/mestelzang.jl")

R0, V0 = 20., 1.
eps0 = 1e-1
ψ(r::Float64)::Float64   = mestel_psi(r,R0,V0,eps0)
dψ(r::Float64)::Float64  = mestel_dpsi_dr(r,R0,V0,eps0)
d2ψ(r::Float64)::Float64 = mestel_ddpsi_ddr(r,R0,V0,eps0)
d3ψ(r::Float64)::Float64 = mestel_d3psi_dr3(r,R0,V0,eps0)
d4ψ(r::Float64)::Float64 = mestel_d4psi_dr4(r,R0,V0,eps0)
# d4ψ(r::Float64)::Float64 = 0.0*r

TOLECC = 0.01
include("../src/Utils/OrbitDefinitions.jl")

if version == "old"
    # Old version
    include("../src/Utils/ComputeEL.jl")
    dEL(a::Float64,e::Float64) = dEdL_from_ae_pot(ψ,dψ,d2ψ,a,e;TOLECC=TOLECC)
else
    # New version 
    include("../src/Utils/ComputeEL_mr.jl")
    dEL(a::Float64,e::Float64) = dEL_from_ae_pot(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;TOLECC=TOLECC)
end

# Array of results 
taba = [0.01, 1.3, 5.4]

emin, emax, δe = 0., 1., 1e-5
tabe = collect(emin:δe:emax)

na = length(taba)
ne = length(tabe)

tabE, tabL = zeros(Float64,na,ne), zeros(Float64,na,ne)
tabdEda, tabdEde = zeros(Float64,na,ne), zeros(Float64,na,ne)
tabdLda, tabdLde = zeros(Float64,na,ne), zeros(Float64,na,ne)

for j = 1:ne
    for i = 1:na
        E, L, ∂E∂a, ∂E∂e, ∂L∂a, ∂L∂e  = dEL(taba[i],tabe[j])

        tabE[i,j], tabL[i,j]        = E, L
        tabdEda[i,j], tabdEde[i,j]  = ∂E∂a, ∂E∂e
        tabdLda[i,j], tabdLde[i,j]  = ∂L∂a, ∂L∂e
    end
end

filename = "ComputeEL_"*version*".hf5"
file = h5open(filename,"w")
write(file,"taba",taba)
write(file,"tabe",tabe)
write(file,"tabE",tabE)
write(file,"tabL",tabL)
write(file,"tabdEda",tabdEda)
write(file,"tabdEde",tabdEde)
write(file,"tabdLda",tabdLda)
write(file,"tabdLde",tabdLde)
close(file)

