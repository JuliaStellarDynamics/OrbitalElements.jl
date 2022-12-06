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
na, ne = 200, 101
amin, amax = 1.e-4, 80.

#####
# Output filename
#####
filename = outdir*"JL_AE_"*modelname*".hf5"


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
    tabJLfromae = zeros(Float64,nae,2)

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
    nJL = nae
    tabaefromJL = zeros(Float64,nJL,2)

    ####
    # Mapping computation
    ####
    Ω₀ = params.Ω₀

    #####
    # (a,e) -> (α,β)
    #####
    println("(a,e) to (J,L) computations :")
    for ae = 1:nae
        print("\r$ae / $nae")
        a, e = tabae[ae,1], tabae[ae,2]
        J, L = OrbitalElements.ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,params)
        if isnan(J) 
            tabJLfromae[ae,1], tabJLfromae[ae,2] = -1., -1.
        elseif isnan(L)
            tabJLfromae[ae,1], tabJLfromae[ae,2] = -2., -1.
        else 
            tabJLfromae[ae,1], tabJLfromae[ae,2] = J, L
        end   
    end
    print("\rDONE                   \n")

    #####
    # (α,β) -> (a,e)
    #####
    println("(J,L) to (a,e) computations :")
    for j = 1:nJL
        print("\r$j / $nJL")
        J, L = tabJLfromae[j,1], tabJLfromae[j,2]
        if L == -1.
            tabaefromJL[j,1], tabaefromJL[j,2] = -1., -1.
        else
            tabaefromJL[j,1], tabaefromJL[j,2] = OrbitalElements.ComputeAEFromActions(ψ,dψ,d2ψ,d3ψ,d4ψ,J,L,params)
        end
    end
    print("\rDONE                   \n")

    
    

    #####
    # Dumping results
    ######
    file = h5open(filename,"w")

    # (a,e) -> (α,β) 
    write(file,"tabAE",tabae)
    write(file,"tabJLfromAE",tabJLfromae)

    # (α,β) -> (a,e)
    write(file,"tabAEfromJL",tabaefromJL)

    close(file)
end


run_test!(filename,
            amin,amax,
            na,ne,
            ψ,dψ,d2ψ,d3ψ,d4ψ,
            params)


# a, e = 10., 1.
# println("a = ",a," ; e = ",e)

# Ω1, Ω2 = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
# println("Ω1 = ",Ω1," ; Ω2 = ",Ω2)

# α, β = Ω1/Ω₀, Ω2/Ω1
# println("α = ",α," ; β = ",β)

# Ω1n, Ω2n = params.Ω₀*α, params.Ω₀*α*β
# println("Ω1n = ",Ω1n," ; Ω2n = ",Ω2n)

# an, en, iter, tol = OrbitalElements.AEFromΩ1Ω2Brute(Ω1n,Ω2n,ψ,dψ,d2ψ,d3ψ,d4ψ,params)
# println("an = ",an," ; en = ",en)
# println("Iterations : ",iter," ; Tolerance : ",tol)

# Ω1l, Ω2l = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,an,en,params)
# println("Ω1l = ",Ω1l," ; Ω2l = ",Ω2l)

# αl, βl = Ω1l/Ω₀, Ω2l/Ω1l
# println("αl = ",αl," ; βl = ",βl)