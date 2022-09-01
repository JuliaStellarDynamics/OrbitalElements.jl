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
# d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψMestel(r,R0,V0,eps0)
# Ω0 = OrbitalElements.Ω0Mestel(R0,V0,eps0)

####
# Isochrone
####
modelname = "Isochrone_mr"
bc, M, G = 1.,1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Ω0 = OrbitalElements.Omega0Isochrone(bc,M,G)

####
# Plummer
####
# modelname = "Plummer_mr"
# bc, M, G = 1.,1.,1.
# ψ(r::Float64)::Float64   = OrbitalElements.plummer_psi(r,bc,M,G)
# dψ(r::Float64)::Float64  = OrbitalElements.plummer_dpsi_dr(r,bc,M,G)
# d2ψ(r::Float64)::Float64 = OrbitalElements.plummer_ddpsi_ddr(r,bc,M,G)
# d3ψ(r::Float64)::Float64 = OrbitalElements.plummer_dddpsi_dddr(r,bc,M,G)
# d4ψ(r::Float64)::Float64 = OrbitalElements.plummer_ddddpsi_ddddr(r,bc,M,G)
# Ω0 = OrbitalElements.plummer_Omega0(bc,M,G)


#####
# Bisections research parameters
#####
rmin = 1.  # For βc research
rmax = 5. # For WminWmax, VminVmax and βc research
#####
# Resonances number (n1,n2) to consider
#####
n1min, n1max = -5, 5
n2min, n2max = -5, 5
#####
# Number of points 
#####
nu, nv = 100, 50
nα, nβ = 100, 50
#####
# Output filename
#####
filename = outdir*"UV_AlphaBeta_"*modelname*".hf5"


function run_test!(filename::String,
                    n1min::Int64,n1max::Int64,n2min::Int64,n2max::Int64,
                    nu::Int64,nv::Int64,nα::Int64,nβ::Int64,
                    ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function;
                    Ω0::Float64=1.,
                    rmin::Float64=1.e-8,rmax::Float64=10000.)

    # β = Ω2 / Ω1 circular function
    βc(α::Float64)::Float64 = OrbitalElements.βcirc(α,dψ,d2ψ,Ω0,rmin=rmin,rmax=rmax)

    nres = (n1max-n1min+1)*(n2max-n2min+1) - 1 # Not adding (0,0) resonance
    tabresnum = zeros(Int64,nres,2)
    tabωminmax = zeros(Float64,nres,2)

    tmp1=1
    for n2 = n2min:n2max
        for n1 = n1min:n1max
            if (n1==0) && (n2==0)
                continue
            end
            tabresnum[tmp1,1], tabresnum[tmp1,2] = n1, n2
            tmp1 += 1
        end
    end

    #####
    # (u,v) points to map to (α,β) -- v points depends on (n1,n2) -> values assigned during computation
    #####
    nuv = nu*nv
    tabu = collect(LinRange(-0.9999,0.9999,nu))
    tabuv = zeros(Float64,nuv,2)
    tabαβfromuv = zeros(Float64,nuv,2)

    tabvminmax = zeros(Float64,nres,nu,2)

    #####
    # (α,β) points to map to (u,v) -- independant of (n1,n2) -> directly filled
    #####
    nαβ = nα*nβ

    #####
    # New version
    #####
    αmin, αmax = OrbitalElements.αminmax(dψ,d2ψ,rmin,rmax,Ω0=Ω0)

    tabα = collect(LinRange(αmin, αmax,nα))
    tabαβ = zeros(Float64,nαβ,2)
    tabuvfromαβ = zeros(Float64,nαβ,2)

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


    file = h5open(filename,"w")

    #####
    # Mapping for each resonance
    #####
    for res = 1:nres

        n1, n2 = tabresnum[res,1], tabresnum[res,2]
        println("Resonance : ($n1,$n2)")
        nameres = "_n1_"*string(n1)*"_n2_"*string(n2)

        # ωmin, ωmax = OrbitalElements.FindWminWmax(n1,n2,dψ,d2ψ,rmax,Ω0) # Old version
        #####
        # New version
        #####
        ωmin, ωmax = OrbitalElements.FindWminWmax(n1,n2,dψ,d2ψ,Ω0=Ω0,rmin=rmin,rmax=rmax)

        tabωminmax[res,1], tabωminmax[res,2] = ωmin, ωmax

        #####
        # (α,β) -> (u,v)
        #####
        for αβ = 1:nαβ
            α, β = tabαβ[αβ,1], tabαβ[αβ,2]
            tabuvfromαβ[αβ,1], tabuvfromαβ[αβ,2] = OrbitalElements.UVFromAlphaBeta(α,β,n1,n2,ωmin,ωmax)
        end
        write(file,"uvfromab"*nameres,tabuvfromαβ)

        #####
        # (u,v) -> (α,β)
        #####
        tmp3=1
        for ju = 1:nu
            uval = tabu[ju]

            #####
            # Old version
            #####
            # vbound = OrbitalElements.FindVbound(n1,n2,dψ,d2ψ,rmax,Ω0)
            # vmin, vmax = OrbitalElements.FindVminVmax(uval,ωmin,ωmax,n1,n2,vbound,βc)
            #####
            # New version
            #####
            vmin, vmax = OrbitalElements.FindVminVmax(uval,n1,n2,dψ,d2ψ,ωmin,ωmax,αmin,αmax,βc,Ω0=Ω0,rmin=rmin,rmax=rmax)
            

            tabvminmax[res,ju,1], tabvminmax[res,ju,2] = vmin, vmax
            dv = (vmax - vmin)/(nv - 1)

            for iv = 1:nv
                vval = vmin + dv*(iv - 1)
                tabαβfromuv[tmp3,1], tabαβfromuv[tmp3,2] = OrbitalElements.AlphaBetaFromUV(uval,vval,n1,n2,ωmin,ωmax)
                tmp3 += 1
            end
        end 
        write(file,"abfromuv"*nameres,tabαβfromuv)
    
    end

    # Resonances
    write(file,"tabres",tabresnum)
    write(file,"tabwminmax",tabωminmax)

    # u (resp. α) interval decomposition
    write(file,"tabu",tabu)
    write(file,"tabAlpha",tabα)

    # Maximal and minimal values for v and β
    write(file,"tabvminmax",tabvminmax)
    write(file,"tabBetaminmax",tabβminmax)

    write(file,"tabAlphaBeta",tabαβ)

    close(file)
end


run_test!(filename,
            n1min,n1max,n2min,n2max,
            nu,nv,nα,nβ,
            ψ,dψ,d2ψ,d3ψ,d4ψ;
            Ω0=Ω0,
            rmin=rmin,rmax=rmax)