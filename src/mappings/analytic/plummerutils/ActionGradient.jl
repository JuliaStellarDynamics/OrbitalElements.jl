

function dPsiEffs(s::Float64, L::Float64;bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    E0 = PlummerE0(bc,M,G)
    #L0 = PlummerL0(bc,M,G)

    if (L != 0.0)
        ds = -E0/s^2 - s*L^2/(bc^2*(s^2-1.0)^2)
        dL = L/(bc^2*(s^2-1.0))
        return ds, dL
    else
        ds = -E0/s^2
        dL = 0.0
        return ds, dL
    end
end


function dJacds(u::Float64, sp::Float64, sa::Float64;bc::Float64=1.,M::Float64=1.,G::Float64=1.)
    Ω₀ = Ω₀Plummer(bc,M,G)

    rp, ra = RpRaFromSpSa(sp,sa,bc)
    E, L   = PlummerELFromSpSa(sp,sa, bc=bc ,M=M,G=G)

    if (sp == 1.0) # Take the limit rp -> 1+
        dPHIdsa, dPHIdLa = dPsiEffs(sa,L, bc=bc ,M=M,G=G)
        dspdE = 0.0
        dsadE = 1.0/dPHIdsa
        dspdL = 0.0
        dsadL = -dPHIdLa/dPHIdsa
    else
        dPHIdsp, dPHIdLp = dPsiEffs(sp,L, bc=bc ,M=M,G=G)
        dPHIdsa, dPHIdLa = dPsiEffs(sa,L, bc=bc ,M=M,G=G)
        dspdE = 1.0/dPHIdsp
        dsadE = 1.0/dPHIdsa
        dspdL = -dPHIdLp/dPHIdsp
        dsadL = -dPHIdLa/dPHIdsa
    end

    fu = u * (1.5 - 0.5*u^2)

    djacdsp = dThetadsp(u, sp, sa,Ω₀)
    djacdsa = dThetadsa(u, sp, sa,Ω₀)

    dsdsp = (1.0 - fu)/2.0
    dsdsa = (1.0 + fu)/2.0

    djacdE = dsadE*djacdsa + dspdE*djacdsp
    djacdL = dsadL*djacdsa + dspdL*djacdsp

    dsdL = dsadL*dsdsa+dspdL*dsdsp

    return djacdE, djacdL, dsdL
end


function GradJrEL(rp::Float64, ra::Float64, nbu::Int64 = 300;bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    Ω₀ = Ω₀Plummer(bc,M,G)

    E, L     = PlummerELFromRpRa(rp,ra, bc=bc ,M=M,G=G)
    sp, sa   = SpSaFromRpRa(rp,ra,bc)
    sma, ecc = AEFromSpSa(sp,sa)

    dJrdE   = 0.0
    dJrdL   = 0.0
    d2JrdE2 = 0.0
    d2JrdEL = 0.0
    d2JrdL2 = 0.0

    # use a midpoint rule to accumulate the derivatives
    for iu=1:nbu
        u = -1.0 + (2/nbu) * (iu-0.5)
        su = SFromUAE(u,sma,ecc)
        ru = RFromS(su,bc)
        xu = ru/bc
        jac = PlummerTheta(u,sp,sa,Ω₀)

        djacdE, djacdL, dsdL = dJacds(u,sp,sa,bc=bc,M=M,G=G)

        # println(su," ",ru," ",xu," ",jac, " ",djacdE)

        dJrdE += jac
        dJrdL += jac/ru^2
        d2JrdE2 += djacdE
        d2JrdEL += djacdL
        d2JrdL2 += jac/ru^2 + L/ru^4*(djacdL*ru^2 - 2.0*bc^2*jac*su*dsdL)
        # println(d2JrdE2)

    end

    dJrdE *= (2/nbu)*1/pi
    dJrdL *= (2/nbu)*(-L/pi)
    d2JrdE2 *= (2/nbu)*1/pi
    d2JrdEL *= (2/nbu)*1/pi
    d2JrdL2 *= (2/nbu)*(-1/pi)

    return dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2
end

function GradJrLlogIntegral(rp::Float64, ra::Float64,
                               nbv::Int64 = 100, eps::Float64=10^(-5),
                               Lcutoff::Float64=0.0005;
                               bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    Ω₀ = Ω₀Plummer(bc,M,G)

    # temporary (ugly) fix until I find the mistake in the function ...
    E, L = PlummerELFromRpRa(rp,ra, bc=bc ,M=M,G=G)
    if (L < Lcutoff)
        L = Lcutoff
        sp, sa = SpSaFromEL(E,L,bc=bc,M=M,G=G)
        rp, ra = RpRaFromSpSa(sp,sa,bc)
    end
    # end of temporary fix

    if (L >= Lcutoff)
        sp, sa = SpSaFromRpRa(rp,ra,bc)
        sma, ecc = AEFromSpSa(sp,sa)
        uminPlusOne = eps*(sp^2-1.0)^2/(3.0*(sa-1.0))

        vmin = log(eps) + 4.0*log(rp) - (log(3.0) + log(sa-1.0))

        dJrdL = 0.0
        d2JrdL2 = 0.0

        for iv=1:nbv
            v = vmin + (log(2.0)-vmin)/nbv * (iv-0.5)
            u = exp(v)-1.0
            su = SFromUAE(u,sma,ecc)
            ru = RFromS(su,bc)
            xu = ru/bc
            jac = PlummerTheta(u,sp,sa,Ω₀)

            djacdE, djacdL, dsdL = dJacds(u,sp,sa,bc=bc,M=M,G=G)

        # println(su," ",ru," ",xu," ",jac, " ",djacdE)


            dJrdL += jac*exp(v)/ru^2

            d2JrdL2 += jac*exp(v)/ru^2 + L*exp(v)/ru^4*(djacdL*ru^2 - 2.0*bc^2*jac*su*dsdL)


        end

        dJrdL *= (log(2.0)-vmin)/nbv * (-L/pi)
        d2JrdL2 *= (log(2.0)-vmin)/nbv * (-1/pi)

        # Leftover integral
        u = uminPlusOne/2.0
        su = SFromUAE(u,sma,ecc)
        ru = RFromS(su,bc)
        jac = PlummerTheta(u,sp,sa,Ω₀)
        djacdE, djacdL, dsdL = dJacds(u,sp,sa,bc=bc,M=M,G=G)

        dJrdL += (uminPlusOne)*jac/ru^2
        d2JrdL2 += (uminPlusOne)*(jac/ru^2 + L/ru^4*(djacdL*ru^2 - 2.0*bc^2*jac*su*dsdL))


        return dJrdL, d2JrdL2

    #error somewhere in there
    else # near-radial orbits
        # maybe we should see what is the exact limit for L=0
        # compute at cutoff for now

        # SOMETHING WRONG HERE
        println(Lcutoff)



        spcut, sacut = SpSaFromEL(E,Lcutoff,bc=bc,M=M,G=G)
        rpcut, racut = RpRaFromSpSa(spcut,sacut,bc)
        smacut, ecccut = AEFromSpSa(spcut,sacut) # this is Plummer (a,e)

        uminPlusOne = eps*(spcut^2-1.0)^2/(3.0*(sacut-1.0))

        vmin = log(eps) + 4.0*log(rpcut) - (log(3.0) + log(sacut-1.0))
        #println(vmin)
        #println("ici")

        dJrdL = 0.0
        d2JrdL2 = 0.0

        for iv=1:nbv
            v = vmin + (log(2.0)-vmin)/nbv * (iv-0.5)
            u = exp(v)-1.0
            su = SFromUAE(u,smacut,ecccut)
            ru = RFromS(su,bc)
            jac = PlummerTheta(u,spcut,sacut,Ω₀)
            djacdE, djacdL, dsdL = dJacds(u,spcut,sacut,bc=bc,M=M,G=G)
            if (rp != 0.0) # not radial
                dJrdL += (uminPlusOne)*jac/ru^2
                d2JrdL2 += (uminPlusOne)*(jac/ru^2 + Lcutoff/ru^4*(djacdL*ru^2 - 2.0*bc^2*jac*su*dsdL))

            end
        end
        dJrdL *= (log(2.0)-vmin)/nbv * (-Lcutoff/pi)
        d2JrdL2 *= (log(2.0)-vmin)/nbv * (-1/pi)

        # Leftover integral
        u   = uminPlusOne/2.0
        su  = SFromUAE(u,smacut,ecccut)
        ru  = RFromS(su,bc)
        jac = PlummerTheta(u,spcut,sacut,Ω₀)

        djacdE, djacdL, dsdL = dJacds(u,spcut,sacut,bc=bc,M=M,G=G)

        dJrdL   += (uminPlusOne)*(-Lcutoff/pi)*jac/ru^2
        d2JrdL2 += (uminPlusOne)*(-1/pi)*(jac/ru^2 + Lcutoff/ru^4*(djacdL*ru^2 - 2.0*bc^2*jac*su*dsdL))

        #println("nonnon")
        if (L == 0.0)
            dJrdL = -0.5
        end

        return dJrdL, d2JrdL2


    end
end

# we should put a cutoff on E near circular orbits
# those should be well defined at circular orbits (cf. isochrone potential)
# ask Mike and Mathieu for near L=0 (Radial) and near Jr=0 (circular)
function GradJrELWrap(E::Float64, L::Float64,
                      nbu::Int64 = 64, eps::Float64=10^(-5), Lcutoff::Float64=0.0005;
                      bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    Ω₀ = Ω₀Plummer(bc,M,G)
    E0 = PlummerE0(bc,M,G)
    L0 = PlummerL0(bc,M,G)

    #println("(E,L) wrap= ",(E,L))
    sp, sa = SpSaFromEL(E,L,bc=bc,M=M,G=G)
    #println("(sp,sa) wrap= ",(sp,sa))
    rp, ra = RpRaFromSpSa(sp,sa,bc)
    #println("(rp,ra) wrap= ",(rp,ra))

    saCutoff = sp + 0.005
    _, raCutoff = RpRaFromSpSa(sp,saCutoff,bc)
    if (sa < saCutoff) # near circular

        # compute exact circular
        dJrdE_circ   = 1/Ω₀ * sp^(5/2)/sqrt(3.0+sp^2)
        dJrdL_circ   = -sp/sqrt(3.0+sp^2)
        d2JrdE2_circ = -(3.0*sp^(7/2)*(-5.0+35.0*sp^2+16.0*sp^4+2.0*sp^6))/(2.0*Ω₀*E0*(3.0+sp^2)^(7/2))
        d2JrdEL_circ = (15.0*sp^2*(sp^2-1.0))/(2.0*E0*(3.0+sp^2)^(7/2))
        d2JrdL2_circ = -3.0/(2.0*L0) * sp^(1/2)*(5.0+7.0*sp^2+4.0*sp^4)/(3.0+sp^2)^(7/2)

        # compute value at cutoff sa

        if (L >= 0.05)
            #println("oui")
            dJrdECut, dJrdLCut, d2JrdE2Cut, d2JrdELCut, d2JrdL2Cut = GradJrEL(rp,raCutoff,nbu, bc=bc ,M=M,G=G)
        else
            dJrdECut, _, d2JrdE2Cut, d2JrdELCut, _ = GradJrEL(rp,raCutoff,nbu, bc=bc ,M=M,G=G)
            dJrdLCut, d2JrdL2Cut = GradJrLlogIntegral(rp,raCutoff,nbu,eps,Lcutoff, bc=bc ,M=M,G=G)
        end

        # taylor expansion near circular orbit (sp,sp)
        # smooth junction at cutoff sa=sp+0.0001

        dJrdE = dJrdE_circ+(dJrdECut-dJrdE_circ)*(sa-sp)/(saCutoff-sp)
        dJrdL = dJrdL_circ+(dJrdLCut-dJrdL_circ)*(sa-sp)/(saCutoff-sp)
        d2JrdE2 = d2JrdE2_circ+(d2JrdE2Cut-d2JrdE2_circ)*(sa-sp)/(saCutoff-sp)
        d2JrdEL = d2JrdEL_circ+(d2JrdELCut-d2JrdEL_circ)*(sa-sp)/(saCutoff-sp)
        d2JrdL2 = d2JrdL2_circ+(d2JrdL2Cut-d2JrdL2_circ)*(sa-sp)/(saCutoff-sp)

        return dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2

    else


        if (L >= 0.05)

            return GradJrEL(rp,ra,nbu, bc=bc ,M=M,G=G)
        else

            dJrdE, _, d2JrdE2, d2JrdEL, _ = GradJrEL(rp,ra,nbu, bc=bc ,M=M,G=G)
            dJrdL, d2JrdL2 = GradJrLlogIntegral(rp,ra,nbu,eps,Lcutoff, bc=bc ,M=M,G=G)
            return dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2
        end
    end
end
