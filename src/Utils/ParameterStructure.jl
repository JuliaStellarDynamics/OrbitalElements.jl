"""
TODO:
1. what to do with default basis and FHT?

"""

struct OrbitsParameters

    Ω₀::Float64
    rmin::Float64
    rmax::Float64
    αmin::Float64
    αmax::Float64

    EDGE::Float64
    TOLECC::Float64
    TOLA::Float64

    NINT::Int64
    FDIFF::Float64

    da::Float64
    de::Float64

    ITERMAX::Int64
    invε::Float64
end


function OrbitsParametersCreate(dψ::Function,d2ψ::Function,Ω₀::Float64;
                                rmin::Float64=1.e-5,rmax::Float64=1.e5,
                                EDGE::Float64=0.01,TOLECC::Float64=0.001,TOLA::Float64=1.0e-4,
                                NINT::Int64=32,
                                FDIFF::Float64=1.0e-8,da::Float64=1.0e-6,de::Float64=1.0e-6,
                                ITERMAX::Int64=100,invε::Float64=1.0e-10)


    # Frequency truncation parameters
    αmin,αmax = OrbitalElements.αminmax(dψ,d2ψ,rmin,rmax,Ω₀)

    return OrbitsParameters(Ω₀,rmin,rmax,αmin,αmax,
                              EDGE,TOLECC,TOLA,NINT,FDIFF,da,de,ITERMAX,invε)
end