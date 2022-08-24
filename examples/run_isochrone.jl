"""
check a bunch of quantities against the isochrone case, so we can be confident we are doing the numerical work correctly!
"""


import OrbitalElements
using Printf

# define easy potentials to pass to frequency calculators
const bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64    = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64   = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64  = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64  = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64  = OrbitalElements.d4ψIsochrone(r,bc,M,G)

Ω₀ = OrbitalElements.Omega0Isochrone(bc,M,G)

x = -1.0
n1,n2 = 1,2
vc = n1*OrbitalElements.Ω1circular(dψ,d2ψ,x) + n2*OrbitalElements.Ω2circular(dψ,x)

println(vc)

# select an (a,e) value for the orbit
a,e = 0.0005, 0.5

# compute rperi and rapo
rp,ra = OrbitalElements.RpRafromAE(a,e)
println("rp=$rp ra=$ra",)

# test frequency computation
Ω₁e,Ω₂e = OrbitalElements.IsochroneOmega12FromAE(a,e,bc,M,G)
Jrr = OrbitalElements.IsochroneJrRpRa(rp,ra,bc,M,G)
#println("Ω₁r=$Ω₁r,Ω₂r=$Ω₂r")

#rcirc1 = OrbitalElements.Omega1circ_to_radius(Ω₁e,dψ,d2ψ)#,Ziter=32,verbose=false)
#println("Ω₁ Bisect r=$rcirc1")

#rcirc0 = OrbitalElements.Omega2circ_to_radius(Ω₂e,dψ)
#println("Ω₂ Bisect r=$rcirc1")

println("truth Ω₁=$Ω₁e,Ω₂=$Ω₂e")

# make a HIGH RES version of the frequencies
Ω₁r,Ω₂r,Jrr = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e,NINT=512,action=true)
println("oldae O1=$Ω₁r O2=$Ω₂r Jr=$Jrr")

#@time Ω₁c,Ω₂c,Jrc = OrbitalElements.HenonThetaFrequenciesRpRa(ψ,dψ,d2ψ,rp,ra,NINT=64,action=true)
#println("theta O1=$Ω₁c O2=$Ω₂c Jr=$Jrc")

@time Ω₁c2,Ω₂c2,Jrc2 = OrbitalElements.HenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT=256,action=true)
println("theta O1=$Ω₁c2 O2=$Ω₂c2 Jr=$Jrc2")

@time Ec,Lc,dEda,dEde,dLda,dLde = OrbitalElements.dELFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)
Em,Lm = OrbitalElements.IsochroneELFromAE(a,e,bc,M,G)
#Em /= OrbitalElements.isochrone_E0()
println("estimated E=$Ec,L=$Lc")
println("true      E=$Em,L=$Lm")


# plot Omega1,Omega2 vs e for the analytic and expansion cases
# is there a problem in the expansion?
# any discontinuities from switching to beta

#=
alpha,beta = Ω₁c/Ω₀,Ω₂c/Ω₁c
println("alpha=$alpha,beta=$beta")

J_EL_ab = OrbitalElements.IsochroneJacELtoAlphaBeta(alpha,beta,bc,M,G)
println("Jacobian(EL,ab):$J_EL_ab")

J_EL_abT = OrbitalElements.JacELToAlphaBetaAE(a,e,ψ,dψ,d2ψ)
println("TJacobian(EL,ab):$J_EL_abT")
=#


# get the numerical frequency derivatives at this point
#f1c,f2c,df1da,df2da,df1de,df2de = OrbitalElements.ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,a,e)
f1c,f2c,df1da,df2da,df1de,df2de = OrbitalElements.ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,a,e,NINT=512,da=1.e-6,de=1.0e-6)


f1h,f2h,df1dah,df1deh,df2dah,df2deh = OrbitalElements.DHenonThetaFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT=128,EDGE=0.2)


# check isochrone numerical diff for frequencies
da = 1.e-6
de = 1.e-6

if e+de > 1.0
    de *= -1.0
end
f1m,f2m = OrbitalElements.IsochroneOmega12FromAE(a,e,bc,M,G)
f1a,f2a = OrbitalElements.IsochroneOmega12FromAE(a+da,e,bc,M,G)
f1e,f2e = OrbitalElements.IsochroneOmega12FromAE(a,e+de,bc,M,G)
df1da2,df1de2,df2da2,df2de2 = (f1a-f1m)/da,(f1e-f1m)/de,(f2a-f2m)/da,(f2e-f2m)/de

println("Compare derivatives:")
println("NDiff : df1da=$df1da,df2da=$df2da,df1de=$df1de,df2de=$df2de")
println("DTheta: df1da=$df1dah,df2da=$df2dah,df1de=$df1deh,df2de=$df2deh")
println("Truth : df1da=$df1da2,df2da=$df2da2,df1de=$df1de2,df2de=$df2de2")


#=
println("df1drp=$df1drp,df2drp=$df2drp,df1dra=$df1dra,df2dra=$df2dra")


J_EL_ae = abs(dEda*dLde - dEde*dLda)
J_o1o2_ae = abs(df1da*df2de - df1de*df2da)
J_o1o2_rpra = abs(df1drp*df2dra - df1dra*df2drp)

tJ_EL_ab = f1c*Ω₀*J_EL_ae/J_o1o2_ae
tJ_EL_abT = f1c*Ω₀*J_EL_ae/J_o1o2_rpra/(2*a)

println("NJacobian(EL,ab):$tJ_EL_ab")
println("NJacobianT(EL,ab):$tJ_EL_abT")

println("Components,$J_EL_ae/$J_o1o2_ae*$f1c*$Ω₀")

# print the analytic and numerical derivatives. who is failing at the boundaries?
println("Checking O1,O2 derivatives")
# frequency derivatives are more sensitive, compared to E,L. mitigation strategies?
# Expansions for frequencies when approaching small values?
da = 1.e-6
de = 1.e-6

if e+de > 1.0
    de *= -1.0
end

f1m,f2m = OrbitalElements.IsochroneOmega12FromAE(a,e,bc,M,G)
f1a,f2a = OrbitalElements.IsochroneOmega12FromAE(a+da,e,bc,M,G)
f1e,f2e = OrbitalElements.IsochroneOmega12FromAE(a,e+de,bc,M,G)
df1da2,df1de2,df2da2,df2de2 = (f1a-f1m)/da,(f1e-f1m)/de,(f2a-f2m)/da,(f2e-f2m)/de
println("Comp derivatives:df1da=$df1da,df1da2=$df1da2")
println("Comp derivatives:df1de=$df1de,df1de2=$df1de2")
println("Comp derivatives:df2da=$df2da,df2da2=$df2da2")
println("Comp derivatives:df2de=$df2de,df2de2=$df2de2")

println("Checking E,L derivatives")
#da = 1.e-5

Em,Lm = OrbitalElements.IsochroneELFromAE(a,e,bc,M,G)
Ea,La = OrbitalElements.IsochroneELFromAE(a+da,e,bc,M,G)
Ee,Le = OrbitalElements.IsochroneELFromAE(a,e+de,bc,M,G)
dEda2,dEde2,dLda2,dLde2 = (Ea-Em)/da,(Ee-Em)/de,(La-Lm)/da,(Le-Lm)/de
println("Comp derivatives:dEda=$dEda,dEda2=$dEda2")
println("Comp derivatives:dEde=$dEde,dEde2=$dEde2")
println("Comp derivatives:dLda=$dLda,dLda2=$dLda2")
println("Comp derivatives:dLde=$dLde,dLde2=$dLde2")


J_EL_aeA = abs(dEda2*dLde2 - dEde2*dLda2)
J_o1o2_aeA = abs(df1da2*df2de2 - df1de2*df2da2)
tJ_EL_abA = f1m*Ω₀*J_EL_aeA/J_o1o2_aeA
println("AJacobian(EL,ab):$tJ_EL_abA")
# do all Jacobians tend to a value when a->0?

=#




#=

# make an accuracy table as a function of number of steps in frequency calculation
open("NINTconvergence.txt","w") do io
    println(io,"truth",",",Ω₁r,",",Ω₂r,",",Jrr)
    for NINT=1:256
        Ω₁c,Ω₂c,Jrc = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e,NINT=NINT,action=true)
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
            rp,ra = OrbitalElements.RpRafromAE(aval,eval)
            Ω₁r,Ω₂r = OrbitalElements.IsochroneOmega12FromRpRa(rp,ra,bc,M,G)
            Jrr = OrbitalElements.isochrone_jr_rpra(rp,ra,bc,M,G)
            #Ω₁r,Ω₂r = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,aval,eval,NINT=1024)
            NINT=8
            o1diff = 1.0
            while (o1diff > TOL)
                #Ω₁c,Ω₂c,Jrc = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,aval,eval,NINT=NINT,action=true)
                Ω₁c,Ω₂c = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,aval,eval,NINT=NINT,TOLECC=TOLECC)#,action=true)
                #Ω₁c,Ω₂c,Jrc = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e,NINT=NINT,action=true)
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
        rp,ra = OrbitalElements.RpRafromAE(aval,eval)
        Ei,Li = OrbitalElements.isochrone_EL_from_rpra(rp,ra,bc,M,G)
        E,L   = OrbitalElements.EL_from_rpra_pot(ψ,dψ,d2ψ,rp,ra)
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
            rp,ra = OrbitalElements.RpRafromAE(aval,eval)
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
=#
