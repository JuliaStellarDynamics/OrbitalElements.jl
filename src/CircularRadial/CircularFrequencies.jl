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

"""Ω1circular(dψ/dr,d²ψ/dr²,a)
radial frequency for circular orbits, from the epicyclic approximation
a is the semi-major axis (equivalent to r for a circular orbit)
"""
function Ω1circular(dψ::Function,
                    d2ψ::Function,
                    a::Float64)

    return sqrt(d2ψ(a) + 3*dψ(a)/a)

end

"""Ω1circular(dψ/dr,d²ψ/dr²,d³ψ/dr³,d⁴ψ/dr⁴,a,e)
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

"""Ω1circular(dψ/dr,d²ψ/dr²,d³ψ/dr³,a,e)
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

"""Ω2circular(dψ/dr,r)
azimuthal frequency for circular orbits, from the epicyclic approximation
"""
function Ω2circular(dψ::Function,a::Float64)

    return sqrt(dψ(a)/a)
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

"""βcircular(dψ/dr,d²ψ/dr²,d³ψ/dr³,a,e)
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


"""βcirc(αcirc, dψ/dr,d²ψ/dr²[, Omega0, rmax])
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
    rcirc = Omega1circ_to_radius(omega1,dψ,d2ψ;rmin=rmin,rmax=rmax)

    omega2 = Omega2_circular(dψ,rcirc)

    return omega2/omega1

end


########################################################################
#
# Radius as a function of circular frequencies (mapping inversion)
#
########################################################################

"""RcircFromΩ1circ(Ω₁,dψ/dr,d²ψ/dr²[, rmin, rmax])
perform backwards mapping from Omega_1 for a circular orbit to radius

can tune [rmin,rmax] for extra optimisation (but not needed)
WARNING: important assumption Ω1_circular is a decreasing function of radius
"""
function RcircFromΩ1circ(ω::Float64,
                        dψ::Function,d2ψ::Function;
                        rmin::Float64=1.0e-8,rmax::Float64=10000.0)


    if ω  <= 0.
        error("Negative circular frequency Ω1 = ",ω)
    end

    rcirc = try bisection(r -> ω - Ω1circular(dψ,d2ψ,r), rmin, rmax) catch;   -1. end
    if (rcirc == -1.) 
        if (0. < ω  < Ω1circular(dψ,d2ψ,rmax))
            return RcircFromΩ1circ(ω,dψ,d2ψ;rmin=rmax,rmax=10*rmax)
        elseif Ω1circular(dψ,d2ψ,rmin) < ω
            error("Too high frequency Ω1 = ",ω)
        else
            error("Unable to find the associated radius of Ω1 = ",ω)
        end
    end

    return rcirc
end


"""RcircFromΩ2circ(Ω₂,dψ/dr[, rmin, rmax])
perform backwards mapping from Omega_2 for a circular orbit to radius

WARNING: important assumption Ω2_circular is a decreasing function of radius
"""
function RcircFromΩ2circ(ω::Float64,
                        dψ::Function;
                        rmin::Float64=1.0e-8,rmax::Float64=10000.0)

    if ω  <= 0.
        error("Negative circular frequency Ω2 = ",ω)
    end

    rcirc = try bisection(r -> ω - Ω2circular(dψ,r), rmin, rmax) catch;  -1. end
    if (rcirc == -1.) 
        if (0. < ω < Ω2circular(dψ,rmax))
            return RcircFromΩ2circ(ω,dψ;rmin=rmax,rmax=10*rmax)
        elseif Ω2circular(dψ,rmin) < ω
            error("Too high frequency Ω2 = ",ω)
        else
            error("Unable to find the associated radius of Ω2 = ",ω)
        end
    end

    return rcirc
end