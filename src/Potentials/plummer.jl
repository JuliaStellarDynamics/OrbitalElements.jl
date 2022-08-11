#=
The plummer potential definitions

  plummer.jl is an example of how to implement a simple new function and derivatives.
  a frequency scaling that creates unity frequencies at the centre is also required.

=#


"""
the plummer potential
"""
function plummer_psi(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=plummer_psi

    the plummer potential
    =#
    rbc = r^2 + bc^2
    return -astronomicalG*M*(sqrt(rbc))^(-1)
end

"""
the plummer potential derivative
"""
function plummer_dpsi_dr(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=plummer_dpsi_dr

    the plummer potential derivative
    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*r*((rbc)^(-3/2))
end

"""
the plummer potential second derivative
"""
function plummer_ddpsi_ddr(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=plummer_ddpsi_ddr

    the plummer potential second derivative
    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*(bc^2 - 2(r^2))*((rbc)^(-5/2))
end

"""
the plummer potential third derivative
"""
function plummer_dddpsi_dddr(r::Float64,bc::Float64=1.0,M::Float64=1.0,astronomicalG::Float64=1.0)::Float64
    rbc = r^2 + bc^2
    term1 = 15*(r^3)/(rbc^(7/2))
    term2 = 9r/(rbc^(5/2))
    return -astronomicalG*M*(-term1 + term2)
end

"""
the plummer potential fourth derivative
"""
function plummer_ddddpsi_ddddr(r::Float64,
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
function plummer_Omega0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=

    plummer frequency scale
    =#
    return 2*sqrt(astronomicalG*M/bc^3)
end
