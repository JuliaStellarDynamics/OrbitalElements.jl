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

"""Ω1circular(dψ,d2ψ,a)
radial frequency for circular orbits, from the epicyclic approximation
a is the semi-major axis (equivalent to r for a circular orbit)
"""
function Ω1circular(dψ::Function,
                    d2ψ::Function,
                    a::Float64)

    if (a == 0.)
        return 2.0*sqrt(abs(d2ψ(0.)))
    else
        return sqrt(d2ψ(a) + 3*dψ(a)/a)
    end
end

"""Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e)
radial frequency for nearly circular orbits, from Taylor expansion
"""
function Ω1circular(dψ::Function,
                    d2ψ::Function,
                    d3ψ::Function,
                    d4ψ::Function,
                    a::Float64,
                    e::Float64)

    # start with the epicylic approximation
    Ω1c = Ω1circular(dψ,d2ψ,a)

    if (e==0.)
        return Ω1c
    end

    # if e is not exactly zero, proceed to the expansion
    dψa, d2ψa, d3ψa, d4ψa = dψ(a), d2ψ(a), d3ψ(a), d4ψ(a)

    # 2nd order Taylor expansion of Omega_1 w.r.t. eccentricity
    zeroorder   = Ω1c
    secondorder =  (
                    -36 * (dψa)^(2)
                    + (a)^3 * (-a * (d3ψa)^(2)
                                + 3*d2ψa*(4*d3ψa + a*d4ψa))
                    + 3 * a * dψa * (12 * d2ψa
                                + a*(20*d3ψa + 3*a*d4ψa))
                    )
                    /
                    (48 * (a)^(2) * (Ω1c)^(3))

    return zeroorder + secondorder * e^2

end

"""Ω1circular(dψ,d2ψ,d3ψ,a,e)
radial frequency for nearly circular orbits, from Taylor expansion
EXCLUDING fourth derivative
"""
function Ω1circular(dψ::Function,
                    d2ψ::Function,
                    d3ψ::Function,
                    a::Float64,
                    e::Float64;
                    FDIFF::Float64=1.e-8)

    # define a numerical fourth derivative
    d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/FDIFF

    return Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e)
end


########################################################################
#
# Azimuthal frequency Ω2 or Frequency ratio β = Ω2 / Ω1
#
########################################################################

"""Ω2circular(dψ,r)
azimuthal frequency for circular orbits, from the epicyclic approximation

@IMPROVE: Taylor expansion in a -> 0+ (need 2nd and 4th derivative)
"""
function Ω2circular(dψ::Function,a::Float64)

    return sqrt(dψ(a)/a)
end

"""Ω2circular(dψ,d2ψ,r)
azimuthal frequency for circular orbits, from the epicyclic approximation
with value at a = 0.
"""
function Ω2circular(dψ::Function,d2ψ::Function,a::Float64)

    if (a == 0.)
        return sqrt(abs(d2ψ(0.)))
    else
        return sqrt(dψ(a)/a)
    end
end

"""
Coefficients of the second-order expansion of β = Ω2/Ω1 near a circular orbit

"""
function βcircular2ndorderExpansionCoefs(ψ::Function,
                                         dψ::Function,
                                         d2ψ::Function,
                                         d3ψ::Function,
                                         d4ψ::Function,
                                         a::Float64)

    # 2nd order Taylor expansion of L
    L0, L1, L2 = Lcirc2ndorderExpansionCoefs(ψ,dψ,d2ψ,d3ψ,a)

    Ω1c = Ω1circular(dψ,d2ψ,a)
    dψa, d2ψa, d3ψa, d4ψa = dψ(a), d2ψ(a), d3ψ(a), d4ψ(a)
    # 2nd order Taylor expansion of Omega_2/(L * Omega_1)
    βoverL0 = 1 / ((a)^(2) * Ω1c)
    βoverL2 =   (
                396 * (dψa)^(2)
                - 3 * a * dψa * (- 100 * d2ψa
                            + 3*a*(4*d3ψa + a*d4ψa))
                + (a)^(2) * (72 * (d2ψa)^(2)
                            + (a)^(2) * (d3ψa)^(2)
                            - a*d2ψa*(4*d3ψa + 3*a*d4ψa))
                )
                /
                (48 * (a)^(4) * (Ω1c)^(5))

    # WARNING: Assumption L1 = Lfirstorder = 0 and βoverLfirstorder = 0
    return L0 * βoverL0, 0., L0 * βoverL0 + L2 * βoverL2
end

"""
Second-order expansion of β = Ω2/Ω1 near a circular orbit
"""
function βcircular(ψ::Function,
                   dψ::Function,
                   d2ψ::Function,
                   d3ψ::Function,
                   d4ψ::Function,
                   a::Float64,
                   e::Float64)

    # compute the Taylor expansion of L
    zeroorder, firstorder, secondorder = βcircular2ndorderExpansionCoefs(ψ,dψ,d2ψ,d3ψ,d4ψ,a)
    return zeroorder + firstorder * e + secondorder * (e)^(2)
end

"""βcircular(dψ,d2ψ,d3ψ,a,e)
β = Ω2/Ω1 for nearly circular orbits, from Taylor expansion
EXCLUDING fourth derivative
"""
function βcircular(ψ::Function,
                   dψ::Function,
                   d2ψ::Function,
                   d3ψ::Function,
                   a::Float64,
                   e::Float64;
                   FDIFF::Float64=1.e-8)

    # define a numerical fourth derivative
    d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/FDIFF

    return βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)
end


"""βcirc(αcirc, dψ,d2ψ[, Omega0, rmax])
return βc(α), the frequency ratio Ω2/Ω1 as a function of α = Ω1/Ω0 .

@IMPROVE: find Omega0 adaptively
"""
function βcirc(αcirc::Float64,
                dψ::Function,d2ψ::Function,
                Ω0::Float64=1.;
                rmin::Float64=1.0e-8,rmax::Float64=10000.)

    # define the circular frequencies
    Ω1 = Ω0 * αcirc
    #rcirc = Omega1circ_to_radius(omega1,dψ,d2ψ,rmax)
    rcirc = RcircFromΩ1circ(Ω1,dψ,d2ψ;rmin=rmin,rmax=rmax)

    Ω2 = Ω2circular(dψ,rcirc)

    return Ω2/Ω1
end


########################################################################
#
# Radius as a function of circular frequencies (mapping inversion)
#
########################################################################

"""RcircFromΩ1circ(Ω₁,dψ,d2ψ[, rmin, rmax])
perform backwards mapping from Omega_1 for a circular orbit to radius

can tune [rmin,rmax] for extra optimisation (but not needed)
@WARNING: important assumption Ω1circular is a decreasing function of radius
"""
function RcircFromΩ1circ(ω::Float64,
                         dψ::Function,d2ψ::Function;
                         rmin::Float64=1.0e-8,rmax::Float64=10000.0,
                         tolx::Float64=1000.0*eps(Float64),tolf::Float64=1000.0*eps(Float64))


    if ω  <= 0.
        error("Negative circular frequency Ω1 = ",ω)
    elseif ω > Ω1circular(dψ,d2ψ,0.)
        error("Too high circular frequency Ω1 = ",ω," > Ω1max = ",Ω1circular(dψ,d2ψ,0.))
    elseif ω == Ω1circular(dψ,d2ψ,0.)
        return 0.
    end

    rcirc = try bisection(r -> ω - Ω1circular(dψ,d2ψ,r),rmin,rmax,tolx=tolx,tolf=tolf) catch;   -1. end

    if (rcirc == -1.)
        if (ω  < Ω1circular(dψ,d2ψ,rmax))
            return RcircFromΩ1circ(ω,dψ,d2ψ;rmin=rmax,rmax=10*rmax,tolx=tolx,tolf=tolf)
        elseif Ω1circular(dψ,d2ψ,rmin) < ω
            return RcircFromΩ1circ(ω,dψ,d2ψ;rmin=rmin/10,rmax=rmin,tolx=tolx,tolf=tolf)
        else
            error("Unable to find the associated radius of Ω1 = ",ω)
        end
    end

    return rcirc
end


"""RcircFromΩ2circ(Ω₂,dψ,d2ψ[, rmin, rmax])
perform backwards mapping from Omega_2 for a circular orbit to radius

@ASSUMPTIONS:
    - Ω2circular is a decreasing function of radius
    - d2ψ used for value at 0.
"""
function RcircFromΩ2circ(ω::Float64,
                        dψ::Function,
                        d2ψ;
                        rmin::Float64=1.0e-8,rmax::Float64=10000.0)

    if ω  <= 0.
        error("Negative circular frequency Ω1 = ",ω)
    elseif ω > Ω2circular(dψ,d2ψ,0.)
        error("Too high circular frequency Ω1 = ",ω," > Ω1max = ",Ω2circular(dψ,d2ψ,0.))
    elseif ω == Ω2circular(dψ,d2ψ,0.)
        return 0.
    end

    rcirc = try bisection(r -> ω - Ω2circular(dψ,r), rmin, rmax) catch;  -1. end
    if (rcirc == -1.)
        if (0. < ω < Ω2circular(dψ,rmax))
            return RcircFromΩ2circ(ω,dψ,d2ψ;rmin=rmax,rmax=10*rmax)
        elseif Ω2circular(dψ,rmin) < ω
            return RcircFromΩ2circ(ω,dψ,d2ψ;rmin=rmin/10,rmax=rmin)
        else
            error("Unable to find the associated radius of Ω2 = ",ω)
        end
    end

    return rcirc
end
