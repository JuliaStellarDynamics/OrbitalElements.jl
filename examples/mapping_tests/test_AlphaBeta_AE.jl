import OrbitalElements
using HDF5

outdir="./"

####
# Mestel
####
# modelname = "Mestel2"
# R0, V0 = 20., 1.
# eps0 = 1e-2
# ψ(r::Float64)::Float64   = OrbitalElements.ψMestel(r,R0,V0,eps0)
# dψ(r::Float64)::Float64  = OrbitalElements.dψMestel(r,R0,V0,eps0)
# d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψMestel(r,R0,V0,eps0)
# d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψMestel(r,R0,V0,eps0)
# Ω0 = OrbitalElements.Ω0Mestel(R0,V0,eps0)

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
modelname = "Plummer"
bc, M, G = 1.,1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.plummer_psi(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.plummer_dpsi_dr(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.plummer_ddpsi_ddr(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.plummer_dddpsi_dddr(r,bc,M,G)
Ω0 = OrbitalElements.plummer_Omega0(bc,M,G)


#####
# Bisections research parameters
#####
rmin, rmax = 1.e-8, 10000. 
#####
# Number of points 
#####
nα, nβ = 5, 2
na, ne = 200, 100
amin, amax = 1.e-3, 10.

#####
# Output filename
#####
filename = outdir*"AlphaBeta_AE_"*modelname*".hf5"


function run_test!(filename::String,
                    nα::Int64,nβ::Int64,
                    na::Int64,ne::Int64,
                    ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function;
                    Ω0::Float64=1.,
                    amin::Float64=1.e-6,amax::Float64=100.,
                    rmin::Float64=1.e-8,rmax::Float64=10000.)

    # β = Ω2 / Ω1 circular function
    βc(α::Float64)::Float64 = OrbitalElements.βcirc(α,dψ,d2ψ,Ω0;rmin=rmin,rmax=rmax)

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
    nαβ = nα*nβ
    tabα = collect(LinRange(0.01,0.99,nα))
    tabαβ = zeros(Float64,nαβ,2)
    tabaefromαβ = zeros(Float64,nαβ,2)

    tabβminmax = zeros(Float64,nα,2)
    # Filling the values
    tmp2=1
    for i = 1:nα
        βmin = 0.5
        βmax = βc(tabα[i])
        dβ = (βmax - βmin)/(nβ-1)
        for j = 1:nβ 
            tabαβ[tmp2,1], tabαβ[tmp2,2] = tabα[i], βmin + (j-1)*dβ
            tmp2 += 1
        end
        tabβminmax[i,1], tabβminmax[i,2] = βmin, βmax
    end

    ####
    # Mapping computation
    ####

    #####
    # (α,β) -> (a,e)
    #####
    println("(alpha,beta) to (a,e) computations :")
    for αβ = 1:nαβ
        print("\r$αβ / $nαβ")
        α, β = tabαβ[αβ,1], tabαβ[αβ,2]
        Ω1, Ω2 = Ω0*α, α*β
        tabaefromαβ[αβ,1], tabaefromαβ[αβ,2] = OrbitalElements.ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,Ω1,Ω2)
    end
    print("\rDONE                   \n")

    #####
    # (a,e) -> (α,β)
    #####
    println("(a,e) to (alpha,beta) computations :")
    for ae = 1:nae
        print("\r$ae / $nae")
        a, e = tabae[ae,1], tabae[ae,2]
        Ω1, Ω2 = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e)
        if isnan(Ω1) 
            tabαβfromae[ae,1], tabαβfromae[ae,2] = 1., 0.
        elseif (Ω1 == 0.) || isnan(Ω2)
            tabαβfromae[ae,1], tabαβfromae[ae,2] = 0., 0.
        else 
            tabαβfromae[ae,1], tabαβfromae[ae,2] = Ω1/Ω0, Ω2/Ω1
        end   
    end
    print("\rDONE                   \n")
    

    #####
    # Dumping results
    ######
    file = h5open(filename,"w")

    # (α,β) -> (a,e)
    write(file,"tabAlphaBeta",tabαβ)
    write(file,"tabAEfromAlphaBeta",tabaefromαβ)

    write(file,"tabAlpha",tabα)
    write(file,"tabBetaminmax",tabβminmax)

    # (a,e) -> (α,β) 
    write(file,"tabAE",tabae)
    write(file,"tabAlphaBetafromAE",tabαβfromae)

    close(file)
end


run_test!(filename,
            nα,nβ,
            na,ne,
            ψ,dψ,d2ψ,d3ψ;
            Ω0=Ω0,
            amin=amin,amax=amax,
            rmin=rmin,rmax=rmax)