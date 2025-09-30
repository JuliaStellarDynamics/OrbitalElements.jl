"""
EdgeHandle.jl

Useful functions for handling edges interpolations.

"""

########################################################################
#
# Checking the variables domain
#
########################################################################
"""
    _checkdomain_ae(a, e)

check for positive semimajor axis and e ∈ [0,1]
"""
function _checkdomain_ae(a, e)
    if a < 0 
        throw(DomainError(a, "Negative semimajor axis")) 
    elseif e < 0 || e > 1 
        throw(DomainError(a, "Unphysical eccentricity not in [0,1]")) 
    end
end

"""
    _checkdomain_rpra(rp, ra)

check for positive pericenter and larger apocenter
"""
function _checkdomain_rpra(rp, ra)
    if rp > ra
        throw(DomainError(rp, "Pericentre larger than apocentre"))
    elseif rp < 0
        throw(DomainError(rp, "Negative pericenter"))
    end
end

########################################################################
#
# Eccentricity tolerance as a function of semi-major axis
#
########################################################################
"""
    _eccentricity_tolerance(ã, tole)

tweak the constant eccentricity tolerance `tole` according to the dimensionless 
semi-major axis.

The semi-major axis has to be rescaled by the system's characteristic size.
For low dimensionless semimajor axis `ã`, increased
tolerance on eccentricity `e` to the domain borders for interpolations.

We effectively want to switch at ``ã*e > TOLE``, i.e., ``tole = TOLE / ã``
where `rc` is the characteristic size of the system.
Constraint eccentricity tolerance to be in [tole,1/3]
@IMPORTANT : 1/3 to prevent infinite recursive loop (leading to StackOverflow Errors)
"""
_eccentricity_tolerance(ã, tole) = min(1/3, max(tole, tole / ã))

########################################################################
#
# Treat any function interpolation in edge cases
#
########################################################################

"""
    _interpolate_edges_ae(fun, ã, e, params)

interpolate the function `fun` close to edges of the dimensionless semimajor axis, 
eccentricity domain `(ã,e)`. 

`fun` has to be given as a function of the dimensionless semimajor axis !
@IMPROVE: Should be a function of the dimensionless semimajor axis
"""
function _interpolate_edges_ae(
    fun::F0,
    a::Float64,
    e::Float64,
    params::OrbitalParameters=OrbitalParameters()
) where {F0 <: Function}

    _checkdomain_ae(a,e)
    # @IMPROVE: TOLA should be dimensionless !
    tola, tole = params.TOLA, _eccentricity_tolerance(a/params.rc, params.TOLECC)
    
    if 0 < a < tola
        # 2nd order interpolation
        # between center (a=0) and value at a=tola and a=2*tola
        # Center:
        a0 = 0.
        f0 = fun(a0, e)
        # ã=tolã:
        atola = tola
        ftola = fun(atola, e)
        # ã=2*tola:
        a2tola = 2tola
        f2tola = fun(a2tola, e)

        return _interpolate_order_2(a, a0, f0, atola, ftola, a2tola, f2tola)
    elseif 0 < e < tole
        # 2nd order interpolation
        # between circular (e=0) and value at e=tole and e=2*tole
        # Circular:
        ecirc = 0.
        fc = fun(a, ecirc)
        # e=tole:
        etole = tole
        ftole = fun(a, etole)
        # e=2*tole:
        e2tole = 2tole
        f2tole = fun(a, e2tole)

        return _interpolate_order_2(e, ecirc, fc, etole, ftole, e2tole, f2tole)
    elseif 1 - tole < e < 1
        # 2nd order interpolation
        # between radial (e=1.) and value at e=1-tole and e=1-2*tole
        # Radial:
        erad = 1.
        fr = fun(a,erad)
        # e=1-tole:
        etole = 1-tole
        ftole = fun(a,etole)
        # e=1-2*tole:
        e2tole = 1-2tole
        f2tole = fun(a,e2tole)

        return _interpolate_order_2(e, e2tole, f2tole, etole, ftole, erad, fr)
    end

    return nothing
end