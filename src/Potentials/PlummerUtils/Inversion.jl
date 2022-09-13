

###
# EVERYTHING FOR FREQUENCY INVERSION
###


function AlphaCirc(beta::Float64)

    num = 0.5 * abs(1.0 - beta^2)^(3/4)
    den = 3.0^(3/4) * beta^(5/2)

    return num/den
end


function TestNearCircular(alpha::Float64, alphac::Float64, sc::Float64;
                bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    E0 = PlummerE0(bc,M,G)
    L0 = PlummerL0(bc,M,G)
    Ω₀ = Ω₀Plummer(bc,M,G)

    d2JrdE2_circ = -(3.0*sc^(7/2)*(-5.0+35.0*sc^2+16.0*sc^4+2.0*sc^6))/(1.0*Ω₀*E0*(3.0+sc^2)^(7/2))
    dE = (alpha-alphac)/(Ω₀*alphac^2*d2JrdE2_circ)
    #println(dE)
    return abs(dE/E0) < 10^(-10)
end


function stepEL(E::Float64, L::Float64, alpha::Float64, alpha_guess::Float64,
                beta::Float64, beta_guess::Float64;
                                bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    Ω₀ = Ω₀Plummer(bc,M,G)

    dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = GradJrELWrap(E,L)

    dalphadE = -Ω₀*alpha_guess^2 * d2JrdE2
    dalphadL = -Ω₀*alpha_guess^2 * d2JrdEL
    dbetadE = -d2JrdEL
    dbetadL = -d2JrdL2

    detJac = dalphadE*dbetadL-dbetadE*dalphadL

    #println(detJac)

    dAlpha = alpha-alpha_guess
    dBeta = beta-beta_guess

    stepE = (dbetadL*dAlpha-dalphadL*dBeta)/detJac
    stepL = (-dbetadE*dAlpha+dalphadE*dBeta)/detJac

    return stepE, stepL
end



function nextEL(E::Float64, L::Float64, stepE::Float64, stepL::Float64;
                bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    E0 = PlummerE0(bc,M,G)
    L0 = PlummerL0(bc,M,G)

    if (E+stepE >= 0.0)
        E = E/2.0
        #println("a")
    elseif (E+stepE <= E0)
        E = (E+E0)/2.0
        #println("b")
    else
        E = E + stepE
        #println("c")
    end


    if (L+stepL > Lc(E,E0,L0))
        L = (L+Lc(E,E0,L0))/2.0
        #println("d")
    elseif (L+stepL <= 0.0)
        L = L/2.0
        #println("e")
    else
        L = L + stepL
        #println("f")
    end
    return E, L
end



function ELRadial(alpha::Float64)
    # L=0.0
    # bissection on alphabetaFromELWrap(E,0.0)
    # bisection(fun, xl, xu)
    L = 0.0
    Emin = -1.0 # alpha=1.0
    Emax = -1.0
    while (PlummerAlphaBetaFromEL(Emax,L)[1] > alpha)
        Emax /= 2.0
    end
    #println((Emin,Emax))
    E = bisection(E->PlummerAlphaBetaFromEL(E,L)[1]-alpha, Emin, Emax)

    alphaguess, betaguess = PlummerAlphaBetaFromEL(E,L)

    return E, L, alpha-alphaguess, 0.0, 1
end

function ELCirc(beta::Float64;
                bc::Float64=1.,M::Float64=1.,G::Float64=1.)
    sc = beta*sqrt(3.0)/sqrt(1.0-beta^2)
    rp, ra = RpRaFromSpSa(sc,sc,bc)
    E, L = PlummerELFromRpRa(rp,ra,bc=bc,M=M,G=G)
    return E, L, 0.0, 0.0, 1
end

# (E,L) From (alpha,beta)
function ELArbitrary(alpha::Float64, beta::Float64, alphac::Float64,
                     nbu::Int64 = 300, eps::Float64=4.0*10^(-10), nbStepMax::Int64=10;
                     bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    E0 = PlummerE0(bc,M,G)
    L0 = PlummerL0(bc,M,G)

    # initial guess (E, L)
    sc = beta*sqrt(3.0)/sqrt(abs(1.0-beta^2))
    rc, _ = RpRaFromSpSa(sc,sc,bc)
    Ec, Lc = PlummerELFromRpRa(rc,rc,bc=bc,M=M,G=G)
    E, L = Ec, Lc

    alphaguess, betaguess = alphac, beta

    # Is the orbit near circular ?
    #dAlpha = dAlpha/dE dE
    # => dE = dAlpha / (dAlpha/dE) evaluated at circular orbits with Lc
    # dAlpha/dE = -Omegam alphac^2 dJr2/dE2 : closed expression for circular orbits

    if (TestNearCircular(alpha,alphac,sc))
        return E, L, alpha-alphaguess, 0.0, 0
    else # Begin Newton's method
        # https://en.wikipedia.org/wiki/Newton%27s_method#k_variables,_k_functions

        nbStep = 0

        while (sqrt((alpha - alphaguess)^2 + (betaguess - beta)^2) > eps)

            stepE, stepL = stepEL(E,L,alpha,alphaguess,beta,betaguess,bc=bc,M=M,G=G)
            E, L = nextEL(E,L,stepE,stepL,bc=bc,M=M,G=G)

            alphaguess, betaguess = PlummerAlphaBetaFromEL(E,L,nbu)
            nbStep += 1

            #println((stepE, stepL,E,L))

            if (nbStep > nbStepMax || (abs(stepE)<abs(E0*eps) && abs(stepL)<abs(L0*eps)) )
                break
            end

        end

        return E, L, alpha - alphaguess, beta - betaguess, nbStep
    end
end
