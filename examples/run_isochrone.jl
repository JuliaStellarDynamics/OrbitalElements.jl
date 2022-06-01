"""
check a bunch of quantities against the isochrone case, so we can be confident we are doing the numerical work correctly!
"""


import OrbitalElements
using Printf

# define easy potentials to pass to frequency calculators
bc, M, G = 1.,1. ,1.
ψ       = r->OrbitalElements.isochrone_psi(r,bc,M,G)
dψdr    = r->OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
d²ψdr²  = r->OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
Ω₀      = OrbitalElements.isochrone_Omega0(bc,M,G)


# select an (a,e) value for the orbit
a,e = 0.01, 0.1

# compute rperi and rapo
rp,ra = OrbitalElements.rpra_from_ae(a,e); @printf("rp=%f ra=%f\n", rp,ra)

# test frequency computation
Ω₁r,Ω₂r = OrbitalElements.isochrone_Omega_1_2(rp,ra,bc,M,G)
Jrr = OrbitalElements.isochrone_jr_rpra(rp,ra,bc,M,G)

# make a HIGH RES version of the frequencies
#Ω₁r,Ω₂r = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e,NINT=1024)
Ω₁c,Ω₂c,Jrc = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e,NINT=32,action=true)
@printf("O1=%f O1guess=%f O2=%f O2guess=%f\n", Ω₁r,Ω₁c,Ω₂r,Ω₂c)

# make an accuracy table as a function of number of steps in frequency calculation
open("NINTconvergence.txt","w") do io
    println(io,"truth",",",Ω₁r,",",Ω₂r,",",Jrr)
    for NINT=1:256
        Ω₁c,Ω₂c,Jrc = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e,NINT=NINT,action=true)
        println(io,NINT,",",Ω₁c,",",Ω₂c,",",Jrc)
    end
end

# make (a,e) grid and search for tolerance
logamin = -5.
logda   = 0.03
na      = 250
emin    = 0.
de      = 0.01
ne      = 100
TOL     = 1.e-4
TOLECC  = 0.03
o1diffmat = zeros(Float64,na,ne)
o2diffmat = zeros(Float64,na,ne)
j1diffmat = zeros(Float64,na,ne)

open("NINTarray.txt","w") do io
    for aindx=1:na
        for eindx=1:ne
            aval = 10.0^(logamin+(aindx-1)*logda)
            eval = (eindx-1)*de
            rp,ra = OrbitalElements.rpra_from_ae(aval,eval)
            Ω₁r,Ω₂r = OrbitalElements.isochrone_Omega_1_2(rp,ra,bc,M,G)
            Jrr = OrbitalElements.isochrone_jr_rpra(rp,ra,bc,M,G)
            #Ω₁r,Ω₂r = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,aval,eval,NINT=1024)
            NINT=8
            o1diff = 1.0
            while (o1diff > TOL)
                #Ω₁c,Ω₂c,Jrc = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,aval,eval,NINT=NINT,action=true)
                Ω₁c,Ω₂c = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,aval,eval,NINT=NINT,TOLECC=TOLECC)#,action=true)
                #Ω₁c,Ω₂c,Jrc = OrbitalElements.compute_frequencies_ae(ψ,dψdr,d²ψdr²,a,e,NINT=NINT,action=true)
                NINT += 1
                o1diff = abs(Ω₁c-Ω₁r)
                o2diff = abs(Ω₂c-Ω₂r)
                j1diff = abs(Jrc-Jrr)
                if NINT>256
                    break
                end
            end
            o1diffmat[aindx,eindx] = NINT-1
            o2diffmat[aindx,eindx] = NINT-1
            j1diffmat[aindx,eindx] = NINT-1
            println(io,aval,",",eval,",",NINT-1)
        end
    end
end


# also do an accuracy grid for re-inverting?


# compute (alpha,beta)
alpha,beta = Ω₁r/Ω₀,Ω₂r/Ω₁r

# invert to get E,L: this can still possibly fail for bad values of sqrts.
E,L = OrbitalElements.isochrone_EL_from_alphabeta(alpha,beta,bc,M,G)
rpcomp,racomp = OrbitalElements.isochrone_rpra_fromEL(E,L,bc,M,G)

# compare to original values
@printf("rp=%f rpguess=%f ra=%f raguess=%f\n", rp,rpcomp,ra,racomp)
acomp,ecomp = OrbitalElements.ae_from_rpra(rpcomp,racomp)
@printf("empircal: a=%f aguess=%f e=%f eguess=%f\n", a,acomp,e,ecomp)

acomp,ecomp = OrbitalElements.isochrone_ae_from_omega1omega2(Ω₁c,Ω₂c,bc,M,G)
@printf("analytical: a=%f aguess=%f e=%f eguess=%f\n", a,acomp,e,ecomp)

# more options for checking...
#
#
# make an accuracy table as a function of number of steps in frequency calculation
aval = 1.
#TOLECC = 0.00001
open("ELconvergence.txt","w") do io
    for eindx=1:1000
        eval = (eindx-1)*0.000001
        rp,ra = OrbitalElements.rpra_from_ae(aval,eval)
        Ei,Li = OrbitalElements.isochrone_EL_from_rpra(rp,ra,bc,M,G)
        E,L   = OrbitalElements.EL_from_rpra_pot(ψ,dψdr,d²ψdr²,rp,ra)
        println(io,eval,",",E-Ei,",",L-Li)
    end
end



"""
o1diffmat = zeros(Float64,na,ne)
o2diffmat = zeros(Float64,na,ne)
j1diffmat = zeros(Float64,na,ne)

open("NINTarray.txt","w") do io
    for aindx=1:na
        for eindx=1:ne
            aval = 10.0^(logamin+(aindx-1)*logda)
            eval = (eindx-1)*de
            rp,ra = OrbitalElements.rpra_from_ae(aval,eval)
            Ei,Li = OrbitalElements.isochrone_EL_from_rpra(rp,ra,bc,M,G)
            NINT=8
            o1diff = 1.0
            while (o1diff > TOL)
                E,L   = OrbitalElements.EL_from_rpra_pot(potential,dpotential,ddpotential,rp,ra)
                NINT += 1
                o1diff = abs(Ω₁c-Ω₁r)
                o2diff = abs(Ω₂c-Ω₂r)
                j1diff = abs(Jrc-Jrr)
                if NINT>256
                    break
                end
            end
            o1diffmat[aindx,eindx] = NINT-1
            o2diffmat[aindx,eindx] = NINT-1
            j1diffmat[aindx,eindx] = NINT-1
            println(io,aval,",",eval,",",NINT-1)
        end
    end
end
"""
