#=
The plummer potential definitions

  plummer.jl is an example of how to implement a simple new function and derivatives.
  a frequency scaling that creates unity frequencies at the centre is also required.

=#


"""
the plummer potential
"""
function ψPlummer(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=ψPlummer

    the plummer potential
    =#
    rbc = r^2 + bc^2
    return -astronomicalG*M*(sqrt(rbc))^(-1)
end

"""
the plummer potential derivative
"""
function dψPlummer(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=dψPlummer

    the plummer potential derivative
    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*r*((rbc)^(-3/2))
end

"""
the plummer potential second derivative
"""
function d2ψPlummer(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=d2ψPlummer

    the plummer potential second derivative
    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*(bc^2 - 2(r^2))*((rbc)^(-5/2))
end

"""
the plummer potential third derivative
"""
function d3ψPlummer(r::Float64,bc::Float64=1.0,M::Float64=1.0,astronomicalG::Float64=1.0)::Float64
    rbc = r^2 + bc^2
    term1 = 15*(r^3)/(rbc^(7/2))
    term2 = 9r/(rbc^(5/2))
    return -astronomicalG*M*(-term1 + term2)
end

"""
the plummer potential fourth derivative
"""
function d4ψPlummer(r::Float64,
                               bc::Float64=1.0,
                               M::Float64=1.0,
                               astronomicalG::Float64=1.0)::Float64
    rbc = r^2 + bc^2
    term1 = 105*(r^4)/(rbc^(9/2))
    term2 = 90*(r^2)/(rbc^(7/2))
    term3 = 9/(rbc^(5/2))
    return -astronomicalG*M*(term1 - term2 + term3)
end


"""
the central frequency for the Plummer potential
"""
function Ω₀Plummer(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=

    plummer frequency scale
    =#
    return 2*sqrt(astronomicalG*M/bc^3)
end
