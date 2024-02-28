"""
EdgeHandle.jl

Useful functions for handling edges interpolations.

"""


########################################################################
#
# Eccentricity tolerance as a function of semi-major axis
#
########################################################################
"""
    EccentricityTolerance(a,rc,TOLECC)

tweak eccentricity tolerance to semi-major axis
"""
function EccentricityTolerance(a::Float64,rc::Float64,TOLECC::Float64)::Float64

    # We effectively want to switch at
    # a*e > rc * TOLECC
    # i.e. tole = rc * TOLECC / a
    # where rc is the characteristic size of the system.
    # + constraint eccentricity tolerance to be in [TOLECC,1/3]
    # Important : 1/3 to prevent infinite recursive loop (leading to StackOverflow Errors)
    return min(1/3, max(TOLECC, TOLECC*rc/a))
end

########################################################################
#
# Treat any function interpolation in edge cases
#
########################################################################

function EdgeHandle(fun::F0,
                    a::Float64,e::Float64,
                    params::OrbitalParameters=OrbitalParameters()) where {F0 <: Function}

    tola, tole = params.TOLA, EccentricityTolerance(a,params.rc,params.TOLECC)
    if (0. < a < tola)
        # 2nd order interpolation
        # between center (a=0) and value at a=tola and a=2*tola
        # Center:
        a0 = 0.
        f0 = fun(a0,e)
        # a=tola:
        atola = tola
        ftola = fun(atola,e)
        # a=2*tola:
        a2tola = 2*tola
        f2tola = fun(a2tola,e)

        # Interpolation
        return Interpolation2ndOrder(a,a0,f0,atola,ftola,a2tola,f2tola)
    elseif (0. < e < tole)
        # 2nd order interpolation
        # between circular (e=0) and value at e=tole and e=2*tole
        # Circular:
        ecirc = 0.
        fc = fun(a,ecirc)
        # e=tole:
        etole = tole
        ftole = fun(a,etole)
        # e=2*tole:
        e2tole = 2*tole
        f2tole = fun(a,e2tole)

        # Interpolation
        return Interpolation2ndOrder(e,ecirc,fc,etole,ftole,e2tole,f2tole)
    elseif ((1.0-tole) < e < 1.)
        # 2nd order interpolation
        # between radial (e=1.) and value at e=1-tole and e=1-2*tole
        # Radial:
        erad = 1.
        fr = fun(a,erad)
        # e=1-tole:
        etole = 1.0-tole
        ftole = fun(a,etole)
        # e=1-2*tole:
        e2tole = 1.0-2*tole
        f2tole = fun(a,e2tole)

        # Interpolation
        return Interpolation2ndOrder(e,e2tole,f2tole,etole,ftole,erad,fr)
    end

    return nothing
end
