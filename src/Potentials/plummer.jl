#=
The plummer potential definitions

  plummer.jl is an example of how to implement a simple new function and derivatives.
  a frequency scaling that creates unity frequencies at the centre is also required.

=#


#=

POTENTIAL AND DERIVATIVES

=#
function plummer_psi(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=plummer_psi

    the plummer potential
    =#
    rbc = r^2 + bc^2
    return -astronomicalG*M*(sqrt(rbc))^(-1)
end

function plummer_dpsi_dr(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=plummer_dpsi_dr

    the plummer potential derivative
    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*r*((rbc)^(-3/2))
end

function plummer_ddpsi_ddr(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=plummer_ddpsi_ddr

    the plummer potential second derivative
    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*(bc^2 - 2(r^2))*((rbc)^(-5/2))
end

function plummer_Omega0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=

    plummer frequency scale
    =#
    return 2*sqrt(astronomicalG*M/bc^3)
end
