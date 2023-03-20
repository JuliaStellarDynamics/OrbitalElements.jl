"""
numerical inversion of (Ω₁,Ω₂) -> (a,e)
by brute-forcing the derivative increments dΩ₁/da, dΩ₁/de, deΩ₂/da, dΩ₂/de

VERBOSE rules:
0: no printing
2: failure cases diagnostics
3: iteration counters

"""


"""
    inverse2Dlinear(a,b,c,d,y1,y2)

Inverse 2D linear system Ax = y, return (0., 0.) if non inversible

A = (a  b)
    (c  d)
"""
function inverse2Dlinear(a::Float64,b::Float64,
                         c::Float64,d::Float64,
                         y1::Float64,y2::Float64)::Tuple{Float64,Float64}

    deta = a*d - b*c
    if deta == 0.
        return 0., 0.
    else
        return (d*y1 - b*y2)/deta, (a*y2 - c*y1)/deta
    end
end


"""
    NextGuessAE(acur,ecur,adir,edir)

Determine the next guess point in the Newton-Rahpson algorithm
given the current point (acur, ecur) and the direction (adir, edir),
dealing with boundary crossing
"""
function NextGuessAE(acur::Float64,ecur::Float64,
                     adir::Float64,edir::Float64,
                     params::OrbitalParameters=OrbitalParameters())

    # If current point at the domain edge
    # and direction towards out.
    if ((ecur == 0.) && (edir < 0.)) || ((ecur == 1.) && (edir > 0.))
        edir = 0.
    elseif (acur == 0.) && (adir < 0.)
        adir = 0.
    end

    # If naive end point (acur + adir, ecur + edir) is outside the safe domain (in the expansion zone):
    # Take the point in the middle of the current point and the border in the increment direction
    # (Dividing the increment length until being in the domain can lead to unexpected behaviour :
    # new point too close to the border where inversion can be impossible).
    tola, tole = params.TOLA, EccentricityTolerance(acur,params.TOLA,params.TOLECC)
    asafemin = 0.01
    esafemin = 0.01
    esafemax = 0.99
    if (acur + adir < tola) || (ecur + edir < tole) || (ecur + edir > 1.0-tole) 
        
        # Fraction of the direction to reach the border
        # acur + afrac * adir = 0.
        afrac = (adir >= 0.) ? Inf : -acur/adir
        if edir == 0.
            efrac = Inf
        elseif edir > 0.
            # Boarder to reach: e = 1.
            # ecur + efrac * edir = 1.
            efrac = (1.0-ecur)/edir
        else
            # Boarder to reach: e = 0.
            # ecur + efrac * edir = 0.
            efrac = -ecur/edir
        end

        # Take the minimal fraction (and not more than 1.0 -> very important !)
        dirfrac = min(min(afrac,efrac),1.0)

        anew = acur + 0.5*dirfrac*adir
        enew = ecur + 0.5*dirfrac*edir
    else
        anew, enew = acur + adir, ecur + edir
    end

    return anew, enew
end



"""
    AEFromΩ1Ω2Brute(Ω₁,Ω₂,ψ,dψ,d2ψ,d3ψ,d4ψ,params)

basic Newton-Raphson algorithm to find (a,e) from (Ω₁,Ω₂) brute force derivatives.
"""
function AEFromΩ1Ω2Brute(Ω₁::Float64,Ω₂::Float64,
                         ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                         params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Int64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # get the circular orbit (maximum radius) for a given Ω₁,Ω₂. use the stronger constraint.
    acirc = RcircFromΩ1circ(Ω₁,dψ,d2ψ,params.rmin,params.rmax)

    # then start from ecc=0.5 and take numerical derivatives
    aguess = acirc
    eguess = 0.5

    f1,f2,df1da,df2da,df1de,df2de = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,aguess,eguess,params)

    tol = (Ω₁ - f1)^2 + (Ω₂ - f2)^2
    if (tol < (params.invε)^2)
        return aguess, eguess, iter, tol
    end

    # 2d Newton Raphson inversion and find new increments
    for iter = 1:params.ITERMAX

        increment1, increment2 = inverse2Dlinear(df1da,df1de,df2da,df2de,Ω₁-f1,Ω₂-f2) 

        # If non invertible 
        if (increment1 == 0.) && (increment2 == 0.)
            break
        end

        # Update guesses
        anew, enew = NextGuessAE(aguess,eguess,increment1,increment2,params)
        aguess, eguess = anew, enew

        f1,f2,df1da,df2da,df1de,df2de = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,aguess,eguess,params)

        tol = (Ω₁ - f1)^2 + (Ω₂ - f2)^2
        if (tol < (params.invε)^2)
            return aguess, eguess, iter, tol
        end
    end

    return aguess,eguess,params.ITERMAX,tol
end


"""
    AEFromJLBrute(J,L,ψ,dψ,d2ψ,d3ψ,params)

basic Newton-Raphson algorithm to find (a,e) from (Jᵣ,L) brute force derivatives.
"""
function AEFromJLBrute(J::Float64,L::Float64,
                       ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                       params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Int64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # get the circular orbit (maximum radius) for a given angular momentum.
    acirc = RcircFromL(L,dψ,params.rmin,params.rmax)

    # then start from ecc=0.5 and take numerical derivatives
    aguess = acirc
    eguess = 0.5

    Jguess, Lguess, dJgda, dLgda, dJgde, dLgde = ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,d3ψ,aguess,eguess,params)

    tol = (Jguess - J)^2 + (Lguess - L)^2
    if (tol < (params.invε)^2)
        return aguess, eguess, iter, tol
    end

    # 2d Newton Raphson inversion and find new increments
    for iter = 1:params.ITERMAX

        increment1, increment2 = inverse2Dlinear(dJgda,dJgde,dLgda,dLgde,J-Jguess,L-Lguess)

        # If non inversible 
        if (increment1 == 0.) && (increment2 == 0.)
            break
        end
        
        # Update guesses
        anew, enew = NextGuess(aguess,eguess,increment1,increment2,params)
        aguess, eguess = anew, enew

        Jguess, Lguess, dJgda, dLgda, dJgde, dLgde = ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,d3ψ,aguess,eguess,params)

        tol = (Jguess - J)^2 + (Lguess - L)^2
        if (tol < (params.invε)^2)
            return aguess, eguess, iter, tol
        end
    end

    return aguess,eguess,params.ITERMAX,tol
end