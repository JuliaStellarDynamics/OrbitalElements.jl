

#####################################
# Default parameter values
#####################################
# recommendations for various default parameters:
const DEFAULT_TOLECC = 0.01
const DEFAULT_TOLA   = 0.1
const DEFAULT_EDGE   = 0.01
const DEFAULT_NINT   = 32
const DEFAULT_RC     = 1.
const DEFAULT_RMIN   = 0.
const DEFAULT_RMAX   = Inf
const DEFAULT_TOL    = 1.e-12
const DEFAULT_DA     = 1.e-5
const DEFAULT_DE     = 1.e-5
const DEFAULT_ITERMAX= 100


"""
    Parameter structure for orbital computations

"""
struct OrbitalParameters

    rc::Float64       # radial scale of the system
    rmin::Float64     # the minimum radius considered for the system
    rmax::Float64     # the maximum radius considered for the system

    EDGE::Float64     # the guard against the edges of [-1,1] integration boundaries
    TOLECC::Float64   # the eccentricity tolerance
    # @IMPROVE: should be dimensionless semimajor axis !
    TOLA::Float64     # the semi-major axis tolerance 

    NINT::Int64       # step size number to compute orbit values (steps from pericentre to apocentre)

    da::Float64       # delta size in semi-major axis for numerical derivatives
    de::Float64       # delta size in eccentricity for numerical derivatives

    ITERMAX::Int64    # number of iterations to minimise inversions of frequency and (E,L) or (a,e)
    invε::Float64     # required accuracy on inversions
end

"""
    OrbitalParameters([,rmin,rmax,EDGE,TOLECC,TOLA,NINT,da,de,ITERMAX,invε])

creates an OrbitalParameters structure with the defined values
(which all have a default value)
"""
function OrbitalParameters(;
    rc::Float64=DEFAULT_RC,
    rmin::Float64=DEFAULT_RMIN,
    rmax::Float64=DEFAULT_RMAX,
    EDGE::Float64=DEFAULT_EDGE,
    TOLECC::Float64=DEFAULT_TOLECC,
    TOLA::Float64=DEFAULT_TOLA,
    NINT::Int64=DEFAULT_NINT,
    da::Float64=DEFAULT_DA,
    de::Float64=DEFAULT_DE,
    ITERMAX::Int64=DEFAULT_ITERMAX,
    invε::Float64=DEFAULT_TOL
)
    return OrbitalParameters(
        rc,
        rmin,
        rmax,
        EDGE,
        TOLECC,
        TOLA,
        NINT,
        da,
        de,
        ITERMAX,
        invε
    )
end
