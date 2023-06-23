"""
    Parameter structure for orbital computations

"""

struct OrbitalParameters

    Ω₀::Float64       # characteristic frequency of the system
    rmin::Float64     # the minimum radius considered for the system
    rmax::Float64     # the maximum radius considered for the system

    rc::Float64       # characteristic size of the system: below this, eccentricity tolerance starts to increase

    EDGE::Float64     # the guard against the edges of [-1,1] integration boundaries
    TOLECC::Float64   # the eccentricity tolerance
    TOLA::Float64     # the semi-major axis tolerance

    NINT::Int64       # step size number to compute orbit values (steps from pericentre to apocentre)

    da::Float64       # delta size in semi-major axis for numerical derivatives
    de::Float64       # delta size in eccentricity for numerical derivatives

    ITERMAX::Int64    # number of iterations to minimise inversions of frequency and (E,L) or (a,e)
    invε::Float64     # required accuracy on inversions
end

"""
    OrbitalParameters([,Ω₀,rmin,rmax,EDGE,TOLECC,TOLA,NINT,da,de,ITERMAX,invε])

creates an OrbitalParameters structure with the defined values (which all have a default value)
"""
function OrbitalParameters(;Ω₀::Float64=DEFAULT_Ω0,
                            rmin::Float64=DEFAULT_RMIN,rmax::Float64=DEFAULT_RMAX,rc::Float64=DEFAULT_RC,
                            EDGE::Float64=DEFAULT_EDGE,TOLECC::Float64=DEFAULT_TOLECC,TOLA::Float64=DEFAULT_TOLA,
                            NINT::Int64=DEFAULT_NINT,
                            da::Float64=DEFAULT_DA,de::Float64=DEFAULT_DE,
                            ITERMAX::Int64=DEFAULT_ITERMAX,invε::Float64=DEFAULT_TOL)

    return OrbitalParameters(Ω₀,rmin,rmax,rc,EDGE,TOLECC,TOLA,NINT,da,de,ITERMAX,invε)
end
