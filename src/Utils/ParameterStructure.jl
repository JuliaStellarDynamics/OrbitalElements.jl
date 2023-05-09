"""
    Parameter structure for orbital computations

"""

struct OrbitalParameters

    Ω₀::Float64
    rmin::Float64
    rmax::Float64

    rc::Float64

    EDGE::Float64
    TOLECC::Float64
    TOLA::Float64

    NINT::Int64

    da::Float64
    de::Float64

    ITERMAX::Int64
    invε::Float64
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
