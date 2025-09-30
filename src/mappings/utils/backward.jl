"""
numerical inversion of some constants -> (a,e)
by brute-forcing the derivative increments dΩ₁/da, dΩ₁/de, deΩ₂/da, dΩ₂/de

"""


"""
    _inverse_2d(a, b, c, d, y1, y2)

Inverse the 2d linear system Ax = y, with
A = (a  b)
    (c  d)

Returns (0., 0.) if non invertible
"""
function _inverse_2d(
    a::Float64,
    b::Float64,
    c::Float64,
    d::Float64,
    y1::Float64,
    y2::Float64
)::Tuple{Float64,Float64}

    deta = a*d - b*c
    if deta == 0. || isnan(deta)
        return 0., 0.
    else
        inc1, inc2 = (d * y1 - b * y2) / deta, (a * y2 - c * y1) / deta
        if isnan(inc1) || isinf(inc1) || isnan(inc2) || isinf(inc2)
            return 0., 0.
        else
            return inc1, inc2
        end
    end
end


"""
    _next_guess_ae(acur, ecur, adir, edir[, params])

Determine the next guess point in the Newton-Rahpson algorithm
given the current point `(acur, ecur)` and the direction `(adir, edir)`,
dealing with boundary crossing

@IMPROVE: Should be a function of the dimensionless semimajor axis
"""
function _next_guess_ae(
    acur::Float64,
    ecur::Float64,
    adir::Float64,
    edir::Float64,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}

    # If current point at the domain edge
    # and direction towards out.
    if (ecur == 0. && edir < 0.) || (ecur == 1. && edir > 0.)
        edir = 0.
    elseif acur == 0. && adir < 0.
        adir = 0.
    end

    # If naive end point (acur + adir, ecur + edir) is outside the safe domain (in the expansion zone):
    # Take the point in the middle of the current point and the border in the increment direction
    # (Dividing the increment length until being in the domain can lead to unexpected behaviour :
    # new point too close to the border where inversion can be impossible).
    # @IMPROVE: For now patched with the parameter rc, but the semimajor axes should be 
    # adimensional
    tola, tole = params.TOLA, _eccentricity_tolerance(acur/params.rc, params.TOLECC)
    if acur + adir < tola || ecur + edir < tole || ecur + edir > 1.0-tole

        # Fraction of the direction to reach the border
        # acur + afrac * adir = 0.
        afrac = adir >= 0. ? Inf : -acur/adir
        if edir == 0.
            efrac = Inf
        elseif edir > 0.
            # Boarder to reach: e = 1.
            # ecur + efrac * edir = 1.
            efrac = (1.0 - ecur) / edir
        else
            # Boarder to reach: e = 0.
            # ecur + efrac * edir = 0.
            efrac = - ecur / edir
        end

        # Take the minimal fraction (and not more than 1.0 -> very important !)
        dirfrac = min(min(afrac, efrac), 1.0)

        anew = acur + dirfrac * adir / 2
        enew = ecur + dirfrac * edir / 2
    else
        anew, enew = acur + adir, ecur + edir
    end

    return anew, enew
end

"""
    _newton_raphson_ae(ainit, einit, v1goal, v2goal, mjacobian[, params])

basic Newton-Raphson algorithm to find `(a,e)`` from `(v1,v2)` goal brute force derivatives.

@ASSUMPTION :
    - The jacobian function is given as
        v1, v2, ∂v1∂a, ∂v2∂a, ∂v1∂e, ∂v2∂e = mjacobian(aguess, eguess)
"""
function _newton_raphson_ae(
    ainit::Float64,
    einit::Float64,
    v1goal::Float64,
    v2goal::Float64,
    mjacobian::F0,
    params::OrbitalParameters=OrbitalParameters()
) where {F0 <: Function}

    # then start from initial guess and take numerical derivatives
    aguess, eguess = ainit, einit
    v1, v2, ∂v1∂a, ∂v2∂a, ∂v1∂e, ∂v2∂e = mjacobian(aguess,eguess)

    tol = (v1goal - v1)^2 + (v2goal - v2)^2
    if (tol < (params.invε)^2)
        return aguess, eguess, 0, tol
    end
    # 2d Newton Raphson inversion and find new increments
    for iter in 1:params.ITERMAX

        adir, edir = _inverse_2d(∂v1∂a, ∂v1∂e, ∂v2∂a, ∂v2∂e, v1goal - v1, v2goal - v2)

        # If non invertible
        if adir == 0. && edir == 0.
            return aguess, eguess, iter, tol
        end

        # Update guesses
        anew, enew = _next_guess_ae(aguess, eguess, adir, edir, params)
        aguess, eguess = anew, enew
        v1, v2, ∂v1∂a, ∂v2∂a, ∂v1∂e, ∂v2∂e = mjacobian(aguess, eguess)

        tol = (v1goal - v1)^2 + (v2goal - v2)^2
        if (tol < (params.invε)^2)
            return aguess, eguess, iter, tol
        end
    end

    return aguess, eguess, params.ITERMAX, tol
end