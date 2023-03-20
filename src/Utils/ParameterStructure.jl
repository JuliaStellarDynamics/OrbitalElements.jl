"""
    Parameter structure for orbital computations

"""

struct OrbitalParameters

    Ω₀::Float64
    rmin::Float64
    rmax::Float64

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
    OrbitalParameters([,Ω₀,rmin,rmax,EDGE,TOLECC,TOLA,NINT,FDIFF,da,de,ITERMAX,invε])

creates an OrbitalParameters structure with the defined values (which all have a default value)
"""
function OrbitalParameters(;Ω₀::Float64=1.0,
                                rmin::Float64=1.e-5,rmax::Float64=1.e5,
                                EDGE::Float64=0.01,TOLECC::Float64=0.01,TOLA::Float64=1.0,
                                NINT::Int64=32,
                                da::Float64=1.0e-6,de::Float64=1.0e-6,
                                ITERMAX::Int64=100,invε::Float64=1.0e-10)

    return OrbitalParameters(Ω₀,rmin,rmax,EDGE,TOLECC,TOLA,NINT,da,de,ITERMAX,invε)
end