"""
CircularFrequencies.jl

Special treatment for circular orbits
  for circular orbits, we can use exact relations with ψ derivatives

Strategies:
-for e=0, use exact mappings for Ω1,Ω2.
-for e~0, use Taylor expansions of Ω1,Ω2.

@IMPROVE, add circular energy mapping

"""


########################################################################
#
# Radial frequency Ω1
#
########################################################################

"""
    Ω1circular(dψ,d2ψ,a)

radial frequency for circular orbits, from the epicyclic approximation
a is the semi-major axis (equivalent to r for a circular orbit)
"""
function Ω1circular(dψ::F1,d2ψ::F2,
                    a::Float64)::Float64 where {F1 <: Function, F2 <: Function}

    if (a == 0.)
        return 2.0*sqrt(abs(d2ψ(0.)))
    else
        return sqrt(d2ψ(a) + 3*dψ(a)/a)
    end
end

########################################################################
#
# Azimuthal frequency Ω2 or Frequency ratio β = Ω2 / Ω1
#
########################################################################

"""
    Ω2circular(dψ,d2ψ,a)

azimuthal frequency for circular orbits, from the epicyclic approximation
with value at a = 0.
"""
function Ω2circular(dψ::F1,d2ψ::F2,
                    a::Float64)::Float64 where {F1 <: Function, F2 <: Function}

    if (a == 0.)
        return sqrt(abs(d2ψ(0.)))
    else
        return sqrt(dψ(a)/a)
    end
end


"""
    βcirc(αcirc,dψ,d2ψ,params)

return βc(α), the frequency ratio Ω2/Ω1 as a function of α = Ω1/Ω₀ .
"""
function βcirc(αcirc::Float64,
               dψ::F1,d2ψ::F2,
               params::OrbitalParameters=OrbitalParameters())::Float64 where {F1 <: Function, F2 <: Function}

    # compute the radial frequency for a circular orbit
    Ω1 = params.Ω₀ * αcirc

    # get the radius corresponding to the circular orbit
    rcirc = RcircFromΩ1circ(Ω1,dψ,d2ψ,params.rmin,params.rmax)

    # get the azimuthal frequency for the radius
    Ω2 = Ω2circular(dψ,d2ψ,rcirc)

    # return β
    return Ω2/Ω1
end


########################################################################
#
# Radius as a function of circular frequencies (mapping inversion)
#
########################################################################

"""
    RadiusFromCircularFrequency(ω,Ω(r),rmin,rmax,tolx,tolf)

perform backwards mapping from frequency for a circular orbit to radius

can tune [rmin,rmax] for extra optimisation (but not needed)
@WARNING: important assumption the frequency is a decreasing function of radius
"""
function RadiusFromCircularFrequency(ω::Float64,
                                     Ωfun::OmF,
                                     rmin::Float64,rmax::Float64,
                                     tolx::Float64=1000.0*eps(Float64),tolf::Float64=1000.0*eps(Float64))::Float64 where {OmF <: Function}

    # check that the input frequency is valid
    if ω  <= 0.
        error("OrbitalElements.Circular.RadiusFromCircularFrequency: Negative circular frequency Ω = $ω")
    elseif ω > Ωfun(0.)
        error("OrbitalElements.Circular.RadiusFromCircularFrequency: Too high circular frequency Ω = $ω > Ω1max = $(Ωfun(0.))")
    elseif ω == Ωfun(0.)
        return 0.
    end

    # use bisection to find the circular orbit radius corresponding to given frequency
    rcirc = try bisection(r -> ω - Ωfun(r),rmin,rmax,tolx=tolx,tolf=tolf) catch;   -1. end

    # check if bisection failed: report why
    if (rcirc == -1.)
        if (ω  < Ωfun(rmax))
            return RadiusFromCircularFrequency(ω,Ωfun,rmax,10*rmax,tolx,tolf)
        elseif Ωfun(rmin) < ω
            return RadiusFromCircularFrequency(ω,Ωfun,rmin/10,rmin,tolx,tolf)
        else
            error("OrbitalElements.Circular.RadiusFromCircularFrequency: Unable to find the associated radius of Ω1 = $ω")
        end
    end

    return rcirc
end

"""
    RcircFromΩ1circ(Ω₁,dψ,d2ψ,rmin,rmax,tolx,tolf)

perform backwards mapping from Omega_1 for a circular orbit to radius

can tune [rmin,rmax] for extra optimisation (but not needed)
@WARNING: important assumption Ω1circular is a decreasing function of radius
"""
function RcircFromΩ1circ(ω::Float64,
                         dψ::F1,d2ψ::F2,
                         rmin::Float64,rmax::Float64,
                         tolx::Float64=1000.0*eps(Float64),tolf::Float64=1000.0*eps(Float64))::Float64 where {F1 <: Function, F2 <: Function}

    # check that the input frequency is valid
    Ωfun(r::Float64)::Float64 = Ω1circular(dψ,d2ψ,r)

    return RadiusFromCircularFrequency(ω,Ωfun,rmin,rmax,tolx,tolf)
end


"""
    RcircFromΩ2circ(Ω₂,dψ,d2ψ,rmin,rmax[,tolx,tolf])

perform backwards mapping from Omega_2 for a circular orbit to radius

@ASSUMPTIONS:
    - Ω2circular is a decreasing function of radius
    - d2ψ used for value at 0.
"""
function RcircFromΩ2circ(ω::Float64,
                         dψ::F1,d2ψ::F2,
                         rmin::Float64,rmax::Float64,
                         tolx::Float64=1000.0*eps(Float64),tolf::Float64=1000.0*eps(Float64))::Float64 where {F1 <: Function, F2 <: Function}

    # check that the input frequency is valid
    Ωfun(r::Float64)::Float64 = Ω2circular(dψ,d2ψ,r)

    return RadiusFromCircularFrequency(ω,Ωfun,rmin,rmax,tolx,tolf)
end