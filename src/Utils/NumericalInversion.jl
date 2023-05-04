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
    if (deta == 0.) || isnan(deta)
        return 0., 0.
    else
        inc1, inc2 = (d*y1 - b*y2)/deta, (a*y2 - c*y1)/deta
        if isnan(inc1) || isinf(inc1) || isnan(inc2) || isinf(inc2)
            return 0., 0.
        else
            return inc1, inc2
        end
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
    tola, tole = params.TOLA, EccentricityTolerance(acur,params.rc,params.TOLECC)
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
    AEFromNewtonRaphson(ainit,einit,v1goal,v2goal,mjacobian,params)

basic Newton-Raphson algorithm to find (a,e) from (v1,v2) goal brute force derivatives.
@ASSUMPTION :
    - The jacobian function is given as 
        v1, v2, ∂v1∂a, ∂v2∂a, ∂v1∂e, ∂v2∂e = mjacobian(aguess,eguess)
"""
function AEFromNewtonRaphson(ainit::Float64,einit::Float64,
                             v1goal::Float64,v2goal::Float64,
                             mjacobian::F0,
                             params::OrbitalParameters=OrbitalParameters()) where {F0 <: Function}

    # then start from initial guess and take numerical derivatives
    aguess, eguess = ainit, einit

    v1, v2, ∂v1∂a, ∂v2∂a, ∂v1∂e, ∂v2∂e = mjacobian(aguess,eguess)

    tol = (v1goal - v1)^2 + (v2goal - v2)^2
    if (tol < (params.invε)^2)
        return aguess, eguess, 0, tol
    end

    # 2d Newton Raphson inversion and find new increments
    for iter = 1:params.ITERMAX

        increment1, increment2 = inverse2Dlinear(∂v1∂a,∂v1∂e,∂v2∂a,∂v2∂e,v1goal-v1,v2goal-v2) 

        # If non invertible 
        if (increment1 == 0.) && (increment2 == 0.)
            return aguess, eguess, iter, tol
        end

        # Update guesses
        anew, enew = NextGuessAE(aguess,eguess,increment1,increment2,params)
        aguess, eguess = anew, enew

        v1, v2, ∂v1∂a, ∂v2∂a, ∂v1∂e, ∂v2∂e = mjacobian(aguess,eguess)

        tol = (v1goal - v1)^2 + (v2goal - v2)^2
        if (tol < (params.invε)^2)
            return aguess, eguess, iter, tol
        end
    end

    return aguess,eguess,params.ITERMAX,tol
end



"""
    AEFromΩ1Ω2Brute(Ω₁,Ω₂,ψ,dψ,d2ψ,params)

basic Newton-Raphson algorithm to find (a,e) from (Ω₁,Ω₂) brute force derivatives.
"""
function AEFromΩ1Ω2Brute(Ω₁::Float64,Ω₂::Float64,
                         ψ::F0,dψ::F1,d2ψ::F2,
                         params::OrbitalParameters=OrbitalParameters()) where {F0 <: Function, F1 <: Function, F2 <: Function}

    # get the circular orbit (maximum radius) for a given Ω₁,Ω₂. use the stronger constraint.
    acirc = RcircFromΩ1circ(Ω₁,dψ,d2ψ,params.rmin,params.rmax)

    # then start from ecc=0.5 and take numerical derivatives
    ainit, einit = acirc, 0.5

    mjac(a::Float64,e::Float64) = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,a,e,params)

    return AEFromNewtonRaphson(ainit,einit,Ω₁,Ω₂,mjac,params)
end


"""
    AEFromJLBrute(J,L,ψ,dψ,d2ψ,params)

basic Newton-Raphson algorithm to find (a,e) from (Jᵣ,L) brute force derivatives.
"""
function AEFromJLBrute(J::Float64,L::Float64,
                       ψ::F0,dψ::F1,
                       params::OrbitalParameters=OrbitalParameters()) where {F0 <: Function, F1 <: Function}

    # Initial guess
    einit = 0.5
    ainitJ = try AfixedEFromJ(J,einit,ψ,dψ,params.rmin,params.rmax,params) catch; 0. end
    ainitL = try AfixedEFromL(L,einit,ψ,dψ,params.rmin,params.rmax,params) catch; 0. end
    ainit = 0.5 * (ainitJ + ainitL)

    mjac(a::Float64,e::Float64) = ComputeActionsAEWithDeriv(ψ,dψ,a,e,params)

    return AEFromNewtonRaphson(ainit,einit,J,L,mjac,params)
end

"""
    AEFromELBrute(J,L,ψ,dψ,d2ψ,params)

basic Newton-Raphson algorithm to find (a,e) from (E,L) brute force derivatives.
"""
function AEFromELBrute(E::Float64,L::Float64,
                       ψ::F0,dψ::F1,
                       params::OrbitalParameters=OrbitalParameters()) where {F0 <: Function, F1 <: Function}

    # Initial guess
    einit = 0.5
    ainitE = try AfixedEFromE(E,einit,ψ,dψ,params.rmin,params.rmax,params) catch; 0. end
    ainitL = try AfixedEFromL(L,einit,ψ,dψ,params.rmin,params.rmax,params) catch; 0. end
    ainit = 0.5 * (ainitE + ainitL)

    mjac(a::Float64,e::Float64) = ComputeELAEWithDeriv(ψ,dψ,a,e,params)

    return AEFromNewtonRaphson(ainit,einit,E,L,mjac,params)
end


########################################################################
#
# Semi-major axis guess for energy, angular momentum and radial action
#
########################################################################
# Guess using Lcirc poorly behave for eccentricity close to one.
# Fix the guess with a given eccentricity

"""
    AfixedEFromL(L,e,ψ,dψ,params,tolx,tolf)
    
perform backwards mapping from energy E for a fixed eccentricity orbit to semi-major axis
can tune [rmin,rmax] for extra optimisation (but not needed)
@WARNING: important assumption E is a increasing function of semi-major axis (at fixed eccentricity)
"""
function AfixedEFromE(E::Float64,e::Float64,
                      ψ::F0,dψ::F1,
                      rmin::Float64,rmax::Float64,
                      params::OrbitalParameters=OrbitalParameters(),
                      tolx::Float64=1000.0*eps(Float64),tolf::Float64=1000.0*eps(Float64))::Float64 where {F0 <:Function, F1 <: Function}

    # check that the input energy is valid
    if E < ψ(0.)
        error("OrbitalElements.Utils.AfixedEFromE: Too low energy E = $E")
    elseif E == ψ(0.)
        return 0.
    elseif E > ψ(Inf)
        error("OrbitalElements.Utils.AfixedEFromE: Too high energy E = $E")
    elseif E == ψ(Inf)
        return Inf
    end

    # use bisection to find the circular orbit radius corresponding to given frequency
    aguess = try bisection(a -> E - EFromAE(ψ,dψ,a,e,params),rmin,rmax,tolx=tolx,tolf=tolf) catch;   -1. end

    # check if bisection failed: report why
    if (aguess == -1.)
        if (EFromAE(ψ,dψ,rmax,e,params) < E)
            return AfixedEFromE(E,e,ψ,dψ,rmax,10*rmax,params,tolx,tolf)
        elseif (E < EFromAE(ψ,dψ,rmin,e,params))
            return AfixedEFromE(E,e,ψ,dψ,rmin/10,rmin,params,tolx,tolf)
        else
            error("OrbitalElements.Utils.AfixedEFromE: Unable to find the associated radius of E = $E and e = $e")
        end
    end

    return aguess
end

"""
    AfixedEFromJ(L,e,ψ,dψ,params,tolx,tolf)
    
perform backwards mapping from radial action J for a fixed eccentricity orbit to semi-major axis
can tune [rmin,rmax] for extra optimisation (but not needed)
@WARNING: important assumption L is a increasing function of semi-major axis (at fixed eccentricity)
"""
function AfixedEFromJ(J::Float64,e::Float64,
                      ψ::F0,dψ::F1,
                      rmin::Float64,rmax::Float64,
                      params::OrbitalParameters=OrbitalParameters(),
                      tolx::Float64=1000.0*eps(Float64),tolf::Float64=1000.0*eps(Float64))::Float64 where {F0 <:Function, F1 <: Function}

    # check that the input energy is valid
    if J < 0.
        error("OrbitalElements.Utils.AfixedEFromJ: Too low radial action J = $J")
    elseif J == 0.
        return 0.
    end

    # use bisection to find the circular orbit radius corresponding to given frequency
    aguess = try bisection(a -> J - HenonJFromAE(ψ,dψ,a,e,params),rmin,rmax,tolx=tolx,tolf=tolf) catch;   -1. end

    # check if bisection failed: report why
    if (aguess == -1.)
        if (HenonJFromAE(ψ,dψ,rmax,e,params) < J)
            return AfixedEFromJ(J,e,ψ,dψ,rmax,10*rmax,params,tolx,tolf)
        elseif (J < HenonJFromAE(ψ,dψ,rmin,e,params))
            return AfixedEFromJ(J,e,ψ,dψ,rmin/10,rmin,params,tolx,tolf)
        else
            error("OrbitalElements.Utils.AfixedEFromJ: Unable to find the associated radius of E = $E and e = $e")
        end
    end

    return aguess
end

"""
    AfixedEFromL(L,e,ψ,dψ,params,tolx,tolf)
    
perform backwards mapping from angular momentum L for a fixed eccentricity orbit to semi-major axis
can tune [rmin,rmax] for extra optimisation (but not needed)
@WARNING: important assumption L is a increasing function of semi-major axis (at fixed eccentricity)
"""
function AfixedEFromL(L::Float64,e::Float64,
                      ψ::F0,dψ::F1,
                      rmin::Float64,rmax::Float64,
                      params::OrbitalParameters=OrbitalParameters(),
                      tolx::Float64=1000.0*eps(Float64),tolf::Float64=1000.0*eps(Float64))::Float64 where {F0 <:Function, F1 <: Function}

    # check that the input frequency is valid
    if L < 0.
        error("OrbitalElements.Utils.AfixedEFromL: Negative angular momentum L = $L")
    elseif L == 0.
        return 0.
    end

    # use bisection to find the circular orbit radius corresponding to given frequency
    aguess = try bisection(a -> L - LFromAE(ψ,dψ,a,e,params),rmin,rmax,tolx=tolx,tolf=tolf) catch;   -1. end

    # check if bisection failed: report why
    if (aguess == -1.)
        if (LFromAE(ψ,dψ,rmax,e,params) < L)
            return AfixedEFromL(L,e,ψ,dψ,rmax,10*rmax,params,tolx,tolf)
        elseif (L < LFromAE(ψ,dψ,rmin,e,params))
            return AfixedEFromL(L,e,ψ,dψ,rmin/10,rmin,params,tolx,tolf)
        else
            error("OrbitalElements.Utils.AfixedEFromL: Unable to find the associated radius of L = $L and e = $e")
        end
    end

    return aguess
end