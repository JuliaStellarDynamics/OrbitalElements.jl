#=
The plummer potential definitions

  plummer.jl is an example of how to implement a simple new function and derivatives.
  a frequency scaling that creates unity frequencies at the centre is also required.

=#


"""
the plummer potential
"""
function ψPlummer(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=ψPlummer

    the plummer potential
    =#
    rbc = r^2 + bc^2
    return -astronomicalG*M*(sqrt(rbc))^(-1)
end

"""
the plummer potential derivative
"""
function dψPlummer(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=dψPlummer

    the plummer potential derivative
    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*r*((rbc)^(-3/2))
end

"""
the plummer potential second derivative
"""
function d2ψPlummer(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=d2ψPlummer

    the plummer potential second derivative
    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*(bc^2 - 2(r^2))*((rbc)^(-5/2))
end

"""
the plummer potential third derivative
"""
function d3ψPlummer(r::Float64,bc::Float64=1.0,M::Float64=1.0,astronomicalG::Float64=1.0)::Float64
    rbc = r^2 + bc^2
    term1 = 15*(r^3)/(rbc^(7/2))
    term2 = 9r/(rbc^(5/2))
    return -astronomicalG*M*(-term1 + term2)
end

"""
the plummer potential fourth derivative
"""
function d4ψPlummer(r::Float64,
                               bc::Float64=1.0,
                               M::Float64=1.0,
                               astronomicalG::Float64=1.0)::Float64
    rbc = r^2 + bc^2
    term1 = 105*(r^4)/(rbc^(9/2))
    term2 = 90*(r^2)/(rbc^(7/2))
    term3 = 9/(rbc^(5/2))
    return -astronomicalG*M*(term1 - term2 + term3)
end


"""
the central frequency for the Plummer potential
"""
function Ω₀Plummer(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    return 2*sqrt(astronomicalG*M/bc^3)
end


"""
Plummer energy scale, from Tep+ 2022 (equation E2)
"""
function PlummerE0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    return -astronomicalG*M/bc
end

"""
Plummer action scale, from Tep+ 2022 (equation E2)
"""
function PlummerL0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    return sqrt(astronomicalG*M*bc)
end



"""
the raw Theta function from Tep et al. 2022, equation F9
"""
function Theta(u::Float64, sp::Float64, sa::Float64,Ω₀::Float64)
    A = sp * (u+2.0)*(u-1.0)^2 - sa*(u-2.0)*(u+1.0)^2
    B = sp * (u^3-3.0*u+6.0) - sa*(u^3-3.0*u-6.0)

    return (1.0/Ω₀ * 3.0/(4.0*sqrt(2.0)) * sqrt(sa*sp*(sa+sp))/sqrt(4.0-u^2)
            * A^(1.5)/sqrt(sa*sp*A + B))
end



function dThetadsp(u::Float64, sp::Float64, sa::Float64,Ω₀::Float64)

    num = (3 *(4 *sp^3 *(-1 + u)^2 *(12 - 3 *u^2 + 2 *u^3 + u^4) +
           sa *sp^2 *(108 - 120 *u - 27 *u^2 + 40 *u^3 + 18 *u^4 - 3 *u^6 +
           3 *sp^2 *(-1 + u)^4 *(2 + u)^2) -
           2 *sa^2 *sp *(-6 - 3 *u + u^3) *(6 - 3 *u + u^3 +
           sp^2 *(2 - 3 *u + u^3)) -
           sa^3 *(-2 + u) *(1 + u)^2 *(6 + 3* u - u^3 + sp^2 *(6 - 3 *u + u^3))))

    den = (8 *sqrt(2)* sp *(sa +
           sp) *sqrt(-(((-4 + u^2) *(sa^2 *sp *(-2 + u) *(1 + u)^2 -
           sp *(6 - 3 *u + u^3) +
           sa *(-6 - 3 *u + u^3 - sp^2 *(2 - 3* u + u^3))))/(
           sa *sp *(sa + sp) *(sa *(-2 + u) *(1 + u)^2 -
           sp *(2 - 3 *u + u^3)))))* (sa^2* sp *(2 + 3 *u - u^3) +
           sp *(6 - 3* u + u^3) + sa *(6 + 3*u - u^3 + sp^2 *(2 - 3 *u + u^3))))

    return 1.0/Ω₀ * num/den
end

function dThetadsa(u::Float64, sp::Float64, sa::Float64,Ω₀::Float64)

    num = (3 *(-4 + u^2) *(3 *sa^4 *sp *(-2 + u)^2 *(1 + u)^4 +
           sp^3 *(-1 + u)^2 *(12 - 3 *u^2 + 2 *u^3 + u^4) -
           2 *sa *sp^2 *(-36 + 9 *u^2 - 6 *u^4 + u^6) -
           2 *sa^3 *(-2 + u) *(1 + u)^2 *(12 + 6 *u - 2 *u^3 +
           sp^2 *(6 - 3 *u + u^3)) +
           sa^2 *sp *(108 + 120 *u - 27 *u^2 - 40*u^3 + 18 *u^4 - 3 *u^6 -
           sp^2 *(-12 + 12 *u + 9 *u^2 - 4 *u^3 - 6 *u^4 + u^6))))

    den = (8 *sqrt(2)*
           sa^2 *sp *(sa + sp)^2 *(sa *(-2 + u) *(1 + u)^2 -
           sp *(2 - 3 *u + u^3)) *(-(((-4 + u^2) *(sa^2 *sp *(-2 + u) *(1 + u)^2 -
           sp *(6 - 3* u + u^3) +
           sa *(-6 - 3* u + u^3 - sp^2 *(2 - 3 *u + u^3))))/(
           sa *sp* (sa + sp) *(sa *(-2 + u) *(1 + u)^2 - sp *(2 - 3 *u + u^3)))))^(3/2))

    return 1.0/Ω₀ * num/den
end


function SpSaFromRpRa(rp::Float64,ra::Float64,bc::Float64)
    sp,sa = sqrt(1+(rp/bc)^2),sqrt(1+(ra/bc)^2)
    return sp,sa
end



"""
the wrapped Theta function
"""
function ΘRpRaPlummer(u::Float64, rp::Float64, ra::Float64, bc::Float64, Ω₀::Float64)
    sp,sa = SpSaFromRpRa(rp,ra,bc)

    return Theta(u,sp,sa,Ω₀)

end


"""
the wrapped Theta derivative function
"""
function dΘRpRaPlummer(u::Float64, rp::Float64, ra::Float64, bc::Float64, Ω₀::Float64)
    sp,sa = SpSaFromRpRa(rp,ra,bc)

    return Theta(u,sp,sa,Ω₀),dThetadsp(u,sp,sa,Ω₀),dThetadsa(u,sp,sa,Ω₀)

end


function sFromUAE(u::Float64, sma::Float64, ecc::Float64)
    fu = u * (1.5 - 0.5*u^2)
    return sma * (1.0 + ecc * fu)
end

function RFroms(s::Float64,bc::Float64)
    return bc * sqrt(abs(s^2 - 1.0))
end

function RFromUAE(u::Float64, a::Float64, e::Float64, bc::Float64)
    return RFroms(sFromUAE(u,a,e),bc)
end

function RFromURpRa(u::Float64, rp::Float64, ra::Float64, bc::Float64)
    sp,sa = SpSaFromRpRa(rp,ra,bc)

    # compute a,e for the anomaly: equation F8 of Tep+ (2022)
    a,e = AEfromRpRa(sp,sa)
    return RFroms(sFromUAE(u,a,e),bc)
end

function PlummerELfromSpSa(sp::Float64, sa::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    E0 = PlummerE0(bc,M,astronomicalG)
    L0 = PlummerL0(bc,M,astronomicalG)

    E = E0/sp - E0*(sa^2-1.0)/(sa*sp*(sa+sp))
    L = L0*sqrt(2.0*(sp^2-1.0)*(sa^2-1.0)/(sa*sp*(sa+sp)))

    return E, L
end

function PlummerELfromRpRa(rp::Float64, ra::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    sp,sa = SpSaFromRpRa(rp,ra,bc)

    return PlummerELfromSpSa(sp,sa)

end

"""Vrad()
SQUARED vr, radial velocity for computing action
as a function of (a,e)
"""
function PlummerVradAE(u::Float64,
                       a::Float64,
                       e::Float64,
                       bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    rp,ra = RpRafromAE(a,e)

    Eval,Lval = PlummerELfromRpRa(rp, ra, bc=bc ,M=M,astronomicalG=astronomicalG)

    r = RFromURpRa(u, rp, ra, bc)

    vrSQ = 2(E-ψPlummer(r,bc,M,astronomicalG)) - (L^2)/(r^2)

    return vrSQ

end

"""Vrad()
SQUARED vr, radial velocity for computing action
as a function of (a,e)
"""
function PlummerVradRpRa(u::Float64,
                         rp::Float64,
                         ra::Float64,
                         bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    # something is weird in here with bc scaling?
    Eval,Lval = PlummerELfromRpRa(rp, ra, bc=bc ,M=M,astronomicalG=astronomicalG)

    r = RFromURpRa(u, rp, ra, bc)

    vrSQ = 2(Eval-ψPlummer(r,bc,M,astronomicalG)) - (Lval^2)/(r^2)

    return vrSQ

end



function PlummerOmega12FromRpRa(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.;NINT=32,action::Bool=false)

    # compute the helpful coordinate for Plummer
    sp,sa = SpSaFromRpRa(rp,ra,bc)

    Ω₀ = Ω₀Plummer(bc,M,astronomicalG)
    Eval,Lval = PlummerELfromSpSa(sp, sa, bc=bc ,M=M,astronomicalG=astronomicalG)

    function u3func(u::Float64)
        # push integration forward on three different quantities: Θ(u),Θ(u)/r^2(u),Θ(u)*vr(u)

        th = ΘRpRaPlummer(u, rp, ra, bc, Ω₀)

        return (th,
                th/(RFromURpRa(u, rp, ra, bc)^2),
                th*PlummerVradRpRa(u,rp,ra,bc,M,astronomicalG))

    end

    accum = UnitarySimpsonIntegration(u3func,K_O=NINT)

    #return the values
    Ω1inv     = (1/pi)*accum[1]
    Ω1        = 1/Ω1inv
    Ω2        = Lval*accum[2]*(1/pi)*Ω1
    actionj   = (1/pi)*accum[3]

    # be careful with Ω2 if near radial: use analytic relation
    #if e>(1-TOLECC)
    #    Ω2 = 0.5*Ω1
    #else
    #    Ω2 = Lval*accum[2]*(1/pi)*Ω1
    #end

    if action
        return Ω1,Ω2,actionj
    else
        return Ω1,Ω2
    end
end



function Sc(tE::Float64)
    t1 = 1.0/(6.0*tE)
    t2 = (1.0+54.0*tE^2)/(216.0*tE^3)
    t3 = (1.0/(4.0*tE)) * sqrt(1.0+1.0/(27.0*tE^2))
    return t1 + cbrt(t2+t3) + cbrt(t2-t3)
end

function tEta(s::Float64, tE::Float64)
    return (s^2 - 1.0)*(1.0/s - tE)
end

function Lc(E::Float64,E0::Float64,L0::Float64)
    tE = E/E0
    if (tE == 1.0)
        return 0.0
    else
        return L0 * sqrt(abs(2*tEta(Sc(tE),tE)))
    end

end

function SpSaFromEL(E::Float64, L::Float64;bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # Use bissection
    # We know that 1 < sp < sc and sc < sa < E0/E for bound non-circular non-radial orbits
    # upper bound should be E0/E + 1 in order to get a proper bracket
    # should find some proper cutoff for quasi circular orbits


    E0 = PlummerE0(bc,M,astronomicalG)
    L0 = PlummerL0(bc,M,astronomicalG)


    tE = E/E0
    tL = L/L0

#    println((E,L))

    if (tE >= 0.0) # bound orbits
        sc = Sc(tE)
        if (L >= Lc(E,E0,L0)) # circular orbit (inequality to take care of small Float errors)
            return sc, sc
        elseif (L == 0.0) # radial orbit
            return 1.0, 1.0/tE
        else # arbitrary orbit
            fct = (s -> tE*s^3 - s^2 + (0.5*tL^2-tE)*s + 1.0)
            if (abs(fct(sc)) <= 10^(-10)) # quasi circular
                return sc, sc
            else # not quasi-circular

                # find sp
                if (abs(fct(1.0)) <= 10.0*eps(Float64))
                    sp = 1.0
                else
                    sp = bisection(fct,1.0,sc)
                end

                #find sa
                if (abs(fct(1.0/tE)) <= 10.0*eps(Float64))
                    sa = 1.0/tE
                else
                    sa = bisection(fct,sc,1.0/tE+1)
                end

                return sp, sa
            end
        end
    else # unbound orbits
        if (L == 0.0) # radial orbits
            return 1.0, Inf
        else # arbitrary orbit
            fct = (s -> tE*s^3 - s^2 + (0.5*tL^2-tE)*s + 1.0)
            sp = bisection(fct,1.0,1.0/(0.5*tL^2-tE))
            return sp, Inf

        end
    end


end
