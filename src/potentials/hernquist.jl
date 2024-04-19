
#=
The Hernquist potential definitions

  hernquist.jl is an example of how to implement a simple new function and derivatives.
  a frequency scaling that creates unity frequencies at the centre is also required.

=#


"""
the plummer potential
"""
function ψHernquist(r::Float64,bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    x = r/bc
    return -(G*M/bc) / sqrt(1.0+x^2)
end

"""
the plummer potential derivative
"""
function dψHernquist(r::Float64,bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    x = r/bc
    # Stable version at infinity (not stable in 0.)
    if x > 1.e5
        return (G*M/(bc^2)) / (x^2 * (sqrt(1.0+x^(-2)))^3)
    end

    return (G*M/(bc^2)) * x / (sqrt(1.0+x^2))^3
end

"""
the plummer potential second derivative
"""
function d2ψHernquist(r::Float64,bc::Float64=1.,M::Float64=1.,G::Float64=1.)
    x = r/bc
    # Stable version at infinity (not stable in 0.)
    if x > 1.e5
        return (G*M/(bc^3)) * ( 1.0 / (sqrt(1.0 + x^2))^3 - 3.0 / (x^3 * (sqrt(1.0 + x^(-2)))^5))
    end

    return (G*M/(bc^3)) * (1.0 - 2.0*(x^2)) / (sqrt(1.0 + x^2))^5
end

"""
the central frequency for the Hernquist potential
"""
function Ω₀Hernquist(bc::Float64=1.,M::Float64=1.,G::Float64=1.)

    return 2*sqrt(G*M/bc^3)
end
