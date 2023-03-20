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

"""
    Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e)

radial frequency for nearly circular orbits, from Taylor expansion
"""
function Ω1circular(dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                    a::Float64,e::Float64)::Float64 where {F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # start with the epicylic approximation
    Ω1c = Ω1circular(dψ,d2ψ,a)

    if (e==0.)
        return Ω1c
    end

    # if e is not exactly zero, proceed to the expansion
    dψa, d2ψa, d3ψa, d4ψa = dψ(a), d2ψ(a), d3ψ(a), d4ψ(a)

    # 2nd order Taylor expansion of Omega_1 w.r.t. eccentricity
    zeroorder   = Ω1c
    # WARNING !!
    # Right operators ar the end of the previous line
    # Very important ! (otherwise some operations might not be done)
    secondorder = (-36 * (dψa)^(2) + 
                    (a)^3 * (-a * (d3ψa)^(2) + 
                    3*d2ψa*(4*d3ψa + a*d4ψa)) + 
                    3 * a * dψa * (12 * d2ψa + 
                    a*(20*d3ψa + 3*a*d4ψa))) / 
                    (48 * (a^2) * (Ω1c^3))

    return zeroorder + secondorder * e^2
end


"""
    dΩ1circular(dψ,d2ψ,d3ψ,a)
"""
function dΩ1circular(dψ::F1,d2ψ::F2,d3ψ::F3,
                     a::Float64)::Float64 where {F1 <: Function, F2 <: Function, F3 <: Function}

    Ω1c = Ω1circular(dψ,d2ψ,a)

    return ((d3ψ(a) + (3/a)*d2ψ(a) - (3/(a^2))*dψ(a))/(2Ω1c))
end

"""
    αcircular(dψ,d2ψ,d3ψ,d4ψ,a,e,Ω₀)

radial frequency for circular orbits normalised by Ω₀
"""
function αcircular(dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                   a::Float64,e::Float64,
                   Ω₀::Float64)::Float64 where {F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    return Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e)/Ω₀
end


########################################################################
#
# Azimuthal frequency Ω2 or Frequency ratio β = Ω2 / Ω1
#
########################################################################

"""
    Ω2circular(dψ,r)

azimuthal frequency for circular orbits, from the epicyclic approximation

@IMPROVE: Taylor expansion in a -> 0+ (need 2nd and 4th derivative)
"""
function Ω2circular(dψ::Function,
                    a::Float64)::Float64

    return sqrt(dψ(a)/a)
end

"""
    Ω2circular(dψ,d2ψ,r)

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
    dΩ2circular(dψ,d2ψ,d3ψ,a)
"""
function dΩ2circular(dψ::F1,d2ψ::F2,d3ψ::F3,
                     a::Float64)::Float64 where {F1 <: Function, F2 <: Function, F3 <: Function}

    Ω2c = Ω2circular(dψ,d2ψ,a)

    return ((1/a)*d2ψ(a) - (1/(a^2)))*dψ(a)/(2Ω2c)

end



"""
    βcircular2ndorderExpansionCoefs(ψ,dψ,d2ψ,d3ψ,d4ψ,a)

Coefficients of the second-order expansion of β = Ω2/Ω1 near a circular orbit

"""
function βcircular2ndorderExpansionCoefs(dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                                         a::Float64)::Tuple{Float64,Float64,Float64} where {F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # 2nd order Taylor expansion of L
    L0, _, L2 = Lcirc2ndorderExpansionCoefs(dψ,d3ψ,a)

    Ω1c = Ω1circular(dψ,d2ψ,a)
    dψa, d2ψa, d3ψa, d4ψa = dψ(a), d2ψ(a), d3ψ(a), d4ψ(a)
    # 2nd order Taylor expansion of Omega_2/(L * Omega_1)
    βoverL0 = 1 / ((a)^(2) * Ω1c)
    # WARNING !!
    # Right parenthesis, division and denominator on the same line
    # Very important ! (otherwise division is not done)
    βoverL2 =   (396 * (dψa)^(2)
                - 3 * a * dψa * (- 100 * d2ψa
                            + 3*a*(4*d3ψa + a*d4ψa))
                + (a)^(2) * (72 * (d2ψa)^(2)
                            + (a)^(2) * (d3ψa)^(2)
                            - a*d2ψa*(4*d3ψa + 3*a*d4ψa))
                )/(48 * (a)^(4) * (Ω1c)^(5))

    # WARNING: Assumption L1 = Lfirstorder = 0 and βoverLfirstorder = 0
    return L0 * βoverL0, 0., L0 * βoverL2 + L2 * βoverL0
end

"""
    βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)

Second-order expansion of β = Ω2/Ω1 near a circular orbit
"""
function βcircular(dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                   a::Float64,e::Float64)::Float64 where {F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # compute the Taylor expansion of L
    zeroorder, firstorder, secondorder = βcircular2ndorderExpansionCoefs(dψ,d2ψ,d3ψ,d4ψ,a)
    return zeroorder + firstorder * e + secondorder * (e)^(2)
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
    RcircFromΩ1circ(Ω₁,dψ,d2ψ,rmin,rmax[,tolx,tolf])

perform backwards mapping from Omega_1 for a circular orbit to radius

can tune [rmin,rmax] for extra optimisation (but not needed)
@WARNING: important assumption Ω1circular is a decreasing function of radius
"""
function RcircFromΩ1circ(ω::Float64,
                         dψ::F1,d2ψ::F2,
                         rmin::Float64,rmax::Float64,
                         tolx::Float64=1000.0*eps(Float64),tolf::Float64=1000.0*eps(Float64))::Float64 where {F1 <: Function, F2 <: Function}

    # check that the input frequency is valid
    if ω  <= 0.
        error("OrbitalElements.Circular.RcircFromΩ1circ: Negative circular frequency Ω1 = $ω")
    elseif ω > Ω1circular(dψ,d2ψ,0.)
        error("OrbitalElements.Circular.RcircFromΩ1circ: Too high circular frequency Ω1 = $ω > Ω1max = $(Ω1circular(dψ,d2ψ,0.))")
    elseif ω == Ω1circular(dψ,d2ψ,0.)
        return 0.
    end

    # use bisection to find the circular orbit radius corresponding to given frequency
    rcirc = try bisection(r -> ω - Ω1circular(dψ,d2ψ,r),rmin,rmax,tolx=tolx,tolf=tolf) catch;   -1. end

    # check if bisection failed: report why
    if (rcirc == -1.)
        if (ω  < Ω1circular(dψ,d2ψ,rmax))
            return RcircFromΩ1circ(ω,dψ,d2ψ,rmax,10*rmax,tolx=tolx,tolf=tolf)
        elseif Ω1circular(dψ,d2ψ,rmin) < ω
            return RcircFromΩ1circ(ω,dψ,d2ψ,rmin/10,rmin,tolx=tolx,tolf=tolf)
        else
            error("OrbitalElements.Circular.RcircFromΩ1circ: Unable to find the associated radius of Ω1 = $ω")
        end
    end

    return rcirc
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

    if ω  <= 0.
        error("OrbitalElements.Circular.RcircFromΩ2circ: Negative circular frequency Ω1 = ",ω)
    elseif ω > Ω2circular(dψ,d2ψ,0.)
        error("OrbitalElements.Circular.RcircFromΩ2circ: Too high circular frequency Ω1 = ",ω," > Ω1max = ",Ω2circular(dψ,d2ψ,0.))
    elseif ω == Ω2circular(dψ,d2ψ,0.)
        return 0.
    end

    rcirc = try bisection(r -> ω - Ω2circular(dψ,r),rmin,rmax,tolx=tolx,tolf=tolf) catch;  -1. end
    if (rcirc == -1.)
        if (0. < ω < Ω2circular(dψ,rmax))
            return RcircFromΩ2circ(ω,dψ,d2ψ,rmax,10*rmax,tolx=tolx,tolf=tolf)
        elseif Ω2circular(dψ,rmin) < ω
            return RcircFromΩ2circ(ω,dψ,d2ψ,rmin/10,rmin,tolx=tolx,tolf=tolf)
        else
            error("OrbitalElements.Circular.RcircFromΩ2circ: Unable to find the associated radius of Ω2 = ",ω)
        end
    end

    return rcirc
end
