import OrbitalElements
using HDF5

outdir="./"

####
# Mestel
####
const modelname = "Mestel"
const R0, V0 = 20., 1.
const ψ(r::Float64)::Float64   = OrbitalElements.ψMestel(r,R0,V0)
const dψ(r::Float64)::Float64  = OrbitalElements.dψMestel(r,R0,V0)
const d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψMestel(r,R0,V0)
const d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψMestel(r,R0,V0)
const d4ψ(r::Float64)::Float64 = OrbitalElements.d3ψMestel(r,R0,V0)
const Ω₀ = OrbitalElements.Ω₀Mestel(R0,V0)

####
# Isochrone
####
# modelname = "Isochrone"
# bc, M, G = 1.,1.,1.
# ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
# dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
# d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
# d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
# Ω0 = OrbitalElements.Omega0Isochrone(bc,M,G)

####
# Plummer
####
# modelname = "Plummer"
# bc, M, G = 1.,1.,1.
# ψ(r::Float64)::Float64   = OrbitalElements.plummer_psi(r,bc,M,G)
# dψ(r::Float64)::Float64  = OrbitalElements.plummer_dpsi_dr(r,bc,M,G)
# d2ψ(r::Float64)::Float64 = OrbitalElements.plummer_ddpsi_ddr(r,bc,M,G)
# d3ψ(r::Float64)::Float64 = OrbitalElements.plummer_dddpsi_dddr(r,bc,M,G)
# Ω0 = OrbitalElements.plummer_Omega0(bc,M,G)

#####
# Parameters
#####
# OrbitalElements parameters
const EDGE = 0.01
const TOLECC = 0.001
# Radii for frequency truncations
const rmin = 0.1
const rmax = 1.e4

const params = OrbitalElements.OrbitsParametersCreate(dψ,d2ψ,Ω₀;rmin=rmin,rmax=rmax,EDGE=EDGE,TOLECC=TOLECC)

#####
# Number of points 
#####
na, ne = 200, 100
amin, amax = 1.e-4, 80.

#####
# Output filename
#####
filename = outdir*"AB_AE_"*modelname*".hf5"


function run_test!(filename::String,
                    amin::Float64,amax::Float64,
                    na::Int64,ne::Int64,
                    ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                    params::OrbitalElements.OrbitsParameters)

    

    #####
    # (a,e) points to map to (α,β) 
    #####
    nae = na*ne
    taba = collect(LinRange(amin,amax,na))
    tabe = collect(LinRange(0.,1.,ne))
    tabae = zeros(Float64,nae,2)
    tabαβfromae = zeros(Float64,nae,2)

    # Filling the values
    tmp1 = 1
    for e = 1:ne
        for a = 1:na
            tabae[tmp1,1], tabae[tmp1,2] = taba[a], tabe[e]
            tmp1 += 1
        end
    end

    #####
    # (α,β) points to map to (a,e)
    #####
    nαβ = nae
    tabaefromαβ = zeros(Float64,nαβ,2)

    ####
    # Mapping computation
    ####
    Ω₀ = params.Ω₀

    #####
    # (a,e) -> (α,β)
    #####
    println("(a,e) to (alpha,beta) computations :")
    for ae = 1:nae
        print("\r$ae / $nae")
        a, e = tabae[ae,1], tabae[ae,2]
        Ω1, Ω2 = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
        if isnan(Ω1) 
            tabαβfromae[ae,1], tabαβfromae[ae,2] = -1., -1.
        elseif (Ω1 == 0.) || isnan(Ω2)
            tabαβfromae[ae,1], tabαβfromae[ae,2] = -2., -1.
        else 
            tabαβfromae[ae,1], tabαβfromae[ae,2] = Ω1/Ω₀, Ω2/Ω1
        end   
    end
    print("\rDONE                   \n")

    #####
    # (α,β) -> (a,e)
    #####
    println("(alpha,beta) to (a,e) computations :")
    for αβ = 1:nαβ
        print("\r$αβ / $nαβ")
        α, β = tabαβfromae[αβ,1], tabαβfromae[αβ,2]
        if β == -1.
            tabaefromαβ[αβ,1], tabaefromαβ[αβ,2] = -1., -1.
        else
            Ω1, Ω2 = Ω₀*α, α*β
            tabaefromαβ[αβ,1], tabaefromαβ[αβ,2] = OrbitalElements.ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,d4ψ,Ω1,Ω2,params)
        end
    end
    print("\rDONE                   \n")

    
    

    #####
    # Dumping results
    ######
    file = h5open(filename,"w")

    # (a,e) -> (α,β) 
    write(file,"tabAE",tabae)
    write(file,"tabAlphaBetafromAE",tabαβfromae)

    # (α,β) -> (a,e)
    write(file,"tabAEfromAlphaBeta",tabaefromαβ)

    close(file)
end


run_test!(filename,
            amin,amax,
            na,ne,
            ψ,dψ,d2ψ,d3ψ,d4ψ,
            params)