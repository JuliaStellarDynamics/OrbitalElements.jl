#=
The isochrone potential definitions

  isochrone.jl is unique in that we can define many of the quantities analytically, which makes for a useful testbed for empirical potentials.

=#

#=

POTENTIAL AND DERIVATIVES

=#
function isochrone_psi(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=isochrone_psi

    the isochrone potential

    =#
    rbc = r^2 + bc^2
    return -astronomicalG*M*(bc+sqrt(rbc))^(-1)
end

function isochrone_dpsi_dr(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=isochrone_dpsi_dr

    the isochrone potential derivative

    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*r*(sqrt(rbc)*(sqrt(rbc)+bc)^2)^(-1)
end

function isochrone_ddpsi_ddr(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=isochrone_ddpsi_ddr

    the isochrone potential second derivative

    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*(- r^(2)*((rbc^(3/2))*(sqrt(rbc)+bc)^2)^(-1)
                          - 2*r^(2)*(rbc*(sqrt(rbc)+bc)^3)^(-1)
                          + (sqrt(rbc)*(sqrt(rbc)+bc)^2)^(-1))
end

#=

SCALING FUNCTIONS

=#
function isochrone_Omega0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # isochrone frequency scale, from Fouvry 21 (appendix G)
    return sqrt(astronomicalG*M/bc^3)
end

function isochrone_E0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # isochrone energy scale, from Fouvry 21 (appendix G)
    return -sqrt(astronomicalG*M/bc)
end

function isochrone_L0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # isochrone action scale, from Fouvry 21 (appendix G)
    return sqrt(astronomicalG*M*bc)
end


#=

FREQUENCY COMPUTATION

=#
function isochrone_omega_ae(rp::Float64,ra::Float64,bc::Float64=1.)
    # isochrone \omega (Fouvry 21 eq. G5)
    sp,sa = spsa_from_rpra(rp,ra)
    return (2/(sp+sa))^(3/2)
end

function isochrone_eta_ae(rp::Float64,ra::Float64,bc::Float64=1.)
    # isochrone \eta (Fouvry 21 eq. G7)
    xp = rp/bc
    xa = ra/bc
    sp,sa = spsa_from_rpra(rp,ra)
    return (1/2)*(1+(xp*xa)/((1+sp)*(1+sa)))
end


function isochrone_Omega_1_2(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # wrapper for returning both \Omega_1 and \Omega_2

    Omega0   = isochrone_Omega0(gc,M,astronomicalG)
    omega_ae = isochrone_omega_ae(rp,ra)
    eta_ae   = isochrone_eta_ae(rp,ra)
    return omega_ae*Omega0,omega_ae*eta_ae*Omega0

end


#=

ENERGY AND ANGULAR MOMENTUM

=#
function isochrone_E_from_rpra(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # isochrone analytic energy, Fouvry 21 G9
    E0 = isochrone_E0()
    sp,sa = spsa_from_rpra(rp,ra)
    return E0/(sp+sa)
end


function isochrone_L_from_rpra(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # isochrone analytic energy, Fouvry 21 G9
    xp = rp/bc
    xa = ra/bc
    L0 = isochrone_L0()
    sp,sa = spsa_from_rpra(rp,ra)
    return sqrt(2)*L0*xp*xa/sqrt((1+sp)*(1+sa)*(sp+sa))
end

function isochrone_dthetadu_from_rpra(r::Float64,u::Float64,rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=

     the isochrone analytic Jacobian, Fouvry 21 G10
    =#
    xp = rp/bc
    xa = ra/bc
    xr = r/bc
    sr = sqrt(1+xr^(2))
    Omega0 = isochrone_Omega0(bc,M,astronomicalG)
    Omega1,Omega2 = isochrone_Omega_1_2(rp,ra,bc,M,astronomicalG)
    sp,sa = spsa_from_rpra(rp,ra)
    return (3/sqrt(2))*(Omega1/Omega0)*(xr/sqrt(4-u^(2)))*(sqrt((sr+sp)*(sr+sa)*(sp+sa))/sqrt((xr+xp)*(xr+xa)))
end



function isochrone_beta_c(alpha::Float64)
    #=

    =#
    return 1/(1+alpha^(2/3))
end


function fanomaly(u::Float64)
    # the Henon anomaly
    u*(1.5 - 0.5*u^(2))
end

function dfdu(u::Float64)
    # the derivative of the Henon anomaly
    1.5*(1.0 - u^(2))
end

function isochrone_drduINVvrfromrpra(rp::Float64,ra::Float64,u::Float64)
    #=


    =#
    bISO = 1.
    Omega0 = 1.
    Sigma, Delta = (ra+rp)*0.5, (ra-rp)*0.5 # Used for the mapping from u
    r = Sigma + Delta*fanomaly(u) # Current value of the radius
    xp, xa, xr = rp/bISO, ra/bISO, r/bISO # Rescaled pericentre, apocentre, and radius
    sqxp, sqxa, sqxr = sqrt(1.0+xp^(2)), sqrt(1.0+xa^(2)), sqrt(1.0+xr^(2)) # Pre-computing the values of sqrt(1+xp^2), sqrt(1+xa^2), and sqrt(1+xr^(2))
    #####
    drduINVvr = (3.0/(sqrt(2.0)))/(Omega0)*xr*sqrt(((sqxr+sqxp)*(sqxr+sqxa)*(sqxp+sqxa))/((xr+xp)*(xr+xa)*(4.0-u^(2))))# Analytical expression of (dr/du)(1/vr), that is always well-posed
    return drduINVvr # Output
end
