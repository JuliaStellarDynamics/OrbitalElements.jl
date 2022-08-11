#=
The isochrone potential definitions

  isochrone.jl is unique in that we can define many of the quantities analytically, which makes for a useful testbed for empirical potentials.

  nearly all quantities for the isochrone potential can be computed analytically

=#

"""isochrone_psi(r[,bc,M,G])

the isochrone potential
"""
function isochrone_psi(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    rbc = r^2 + bc^2
    return -astronomicalG*M*(bc+sqrt(rbc))^(-1)
end

"""isochrone_dpsi_dr(r[,bc,M,G])

the isochrone potential derivative
"""
function isochrone_dpsi_dr(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    rbc = r^2 + bc^2
    return astronomicalG*M*r*(sqrt(rbc)*(sqrt(rbc)+bc)^2)^(-1)
end

"""isochrone_ddpsi_ddr(r[,bc,M,G])

the isochrone potential second derivative
"""
function isochrone_ddpsi_ddr(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=
    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*(- r^(2)*((rbc^(3/2))*(sqrt(rbc)+bc)^2)^(-1)
                          - 2*r^(2)*(rbc*(sqrt(rbc)+bc)^3)^(-1)
                          + (sqrt(rbc)*(sqrt(rbc)+bc)^2)^(-1))
end

"""isochrone_dddpsi_dddr(r[,bc,M,G])

the isochrone potential third derivative
"""
function isochrone_dddpsi_dddr(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=
    =#
    r2      = r^2
    r3      = r^3
    rbc     = r2 + bc^2
    sqrbc   = sqrt(rbc)
    bcsqrbc = bc+sqrbc
    return -astronomicalG*M*(-(6r3)/((rbc^(3/2))*(bcsqrbc^4))
                             -(6r3)/((rbc^2)*(bcsqrbc^3))
                             +(6r)/(rbc*(bcsqrbc^3))
                             -(3r3)/((rbc^(5/2))*(bcsqrbc^2))
                             +(3r)/((rbc^(3/2))*bcsqrbc^2))
end

"""isochrone_ddddpsi_ddddr(r[,bc,M,G])

the isochrone potential fourth derivative
"""
function isochrone_ddddpsi_ddddr(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=
    =#
    r2      = r^2
    r3      = r^3
    r4      = r^4
    rbc     = r2 + bc^2
    sqrbc   = sqrt(rbc)
    bcsqrbc = bc+sqrbc
    return -astronomicalG*M*(+(24r4)/((rbc^2)*(bcsqrbc^5))
                             +(36r4)/((rbc^(5/2))*(bcsqrbc^4))
                             -(36r2)/((rbc^(3/2))*(bcsqrbc^4))
                             +(30r4)/((rbc^(3))*(bcsqrbc^3))
                             -(36r2)/((rbc^(2))*(bcsqrbc^3))
                             +(6)/((rbc)*(bcsqrbc^3))
                             +(15r4)/((rbc^(7/2))*(bcsqrbc^2))
                             -(18r2)/((rbc^(5/2))*(bcsqrbc^2))
                             +(3)/((rbc^(3/2))*(bcsqrbc^2)))
end



"""
isochrone frequency scale, from Fouvry 21 (appendix G)
"""
function isochrone_Omega0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    return sqrt(astronomicalG*M/bc^3)
end

"""
isochrone energy scale, from Fouvry 21 (appendix G)
Emin = -GM/(2bc)
"""
function isochrone_E0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #return -sqrt(astronomicalG*M/bc)
    return -astronomicalG*M/bc
end

"""
isochrone action scale, from Fouvry 21 (appendix G)
"""
function isochrone_L0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    return sqrt(astronomicalG*M*bc)
end


"""
isochrone reduced coordinates for pericentre and apocentre
"""
function IsochroneSpSaFromRpRa(rp::Float64,ra::Float64,bc::Float64=1.)
    xp = rp/bc
    xa = ra/bc
   return sqrt(1+xp^2),sqrt(1+xa^2)
end

"""
compute the radial action
(Fouvry 21 eq. G3)
"""
function isochrone_jr_rpra(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    E = isochrone_E_from_rpra(rp,ra,bc,M,astronomicalG)
    L = isochrone_L_from_rpra(rp,ra,bc,M,astronomicalG)
    return (astronomicalG*M/sqrt(-2E)) - 0.5 * (L + sqrt(L*L + 4*astronomicalG*M*bc))
end

"""
compute the dimensionless function for Omega1
(Fouvry 21 eq. G5)
"""
function IsochroneAlphaRpRa(rp::Float64,ra::Float64,bc::Float64=1.)
    sp,sa = IsochroneSpSaFromRpRa(rp,ra,bc)
    return (2/(sp+sa))^(3/2)
end

"""
compute the dimensionless function for Omega1
(Fouvry 21 eq. G5)
"""
function IsochroneAlphaAE(a::Float64,ecc::Float64,bc::Float64=1.)
    sp,sa = IsochroneSpSaFromRpRa(a*(1-ecc),a*(1+ecc),bc)
    return (2/(sp+sa))^(3/2)
end

"""
compute the dimensionless function for Omega2 from (rp,ra)
(Fouvry 21 eq. G7)
"""
function IsochroneBetaRpRa(rp::Float64,ra::Float64,bc::Float64=1.)
    xp = rp/bc
    xa = ra/bc
    sp,sa = IsochroneSpSaFromRpRa(rp,ra,bc)
    return (1/2)*(1+(xp*xa)/((1+sp)*(1+sa)))
end

"""
compute the dimensionless function for Omega2 from (a,e)
(Fouvry 21 eq. G7)
"""
function IsochroneBetaAE(a::Float64,ecc::Float64,bc::Float64=1.)
    xp = (a*(1-ecc))/bc
    xa = (a*(1+ecc))/bc
    sp,sa = IsochroneSpSaFromRpRa(a*(1-ecc),a*(1+ecc),bc)
    return (1/2)*(1+(xp*xa)/((1+sp)*(1+sa)))
end


"""
analytic function to return isochrone frequencies from (rp,ra)
"""
function IsochroneOmega12FromRpRa(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    Omega0   = isochrone_Omega0(bc,M,astronomicalG)
    omega_ae = IsochroneAlphaRpRa(rp,ra,bc)
    eta_ae   = IsochroneBetaRpRa(rp,ra,bc)
    return omega_ae*Omega0,omega_ae*eta_ae*Omega0

end

"""
analytic function to return isochrone frequencies from (a,e)
"""
function IsochroneOmega12FromAE(a::Float64,ecc::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    Omega0   = isochrone_Omega0(bc,M,astronomicalG)
    omega_ae = IsochroneAlphaAE(a,ecc,bc)
    eta_ae   = IsochroneBetaAE(a,ecc,bc)
    return omega_ae*Omega0,omega_ae*eta_ae*Omega0
end

"""
inversion of EL -> alpha,beta function
"""
function isochrone_EL_from_alphabeta(alpha::Float64,beta::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    scaleAction = isochrone_L0(bc,M,astronomicalG)
    E = 0.5*scaleEnergy*(alpha)^(2/3) # Value of the energy
    L = scaleAction*(2.0*beta-1.0)/(sqrt(beta*(1.0-beta))) # Value of the angular momentum
    return E, L # Output
end

"""

@IMPROVE, uses a floor to avoid any sqrt problems with circular orbits
"""
function isochrone_rpra_fromEL(E::Float64,L::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    scaleAction = isochrone_L0(bc,M,astronomicalG)
    xc = ((0.5*scaleEnergy)/E) - 1.0                                 # Value of xc
    ecc_ov = sqrt(max(0,1.0 - (L/scaleAction)^(2)*(1.0/xc)*(1.0+(1.0/xc)))) # Value of overline{e}
    xp = sqrt(max(0,(2.0 + xc*(1.0-ecc_ov))*(xc*(1.0-ecc_ov))))             # Value of xp
    xa = sqrt(max(0,(2.0 + xc*(1.0+ecc_ov))*(xc*(1.0+ecc_ov))))             # Value of xa
    rp, ra = xp*bc, xa*bc                                            # Value of rp,ra
    return rp, ra # Output
end

"""
function to wrap (alpha,beta)->(E,L)->(rp,ra)->(a,e) conversions for isochrone
"""
function isochrone_ae_from_omega1omega2(omega1::Float64,omega2::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    Omega0= isochrone_Omega0(bc,M,astronomicalG)
    E,L   = isochrone_EL_from_alphabeta(omega1/Omega0,omega2/omega1,bc,M,astronomicalG)
    rp,ra = isochrone_rpra_fromEL(E,L,bc,M,astronomicalG)
    a,e   = ae_from_rpra(rp,ra)
    return a,e
end

"""
energy from isochrone model, using rpra
Fouvry 21 G9
"""
function isochrone_E_from_rpra(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    sp,sa       = IsochroneSpSaFromRpRa(rp,ra,bc)
    return scaleEnergy/(sp+sa)
end

"""
angular momentum from isochrone model, using rpra
"""
function isochrone_L_from_rpra(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # isochrone analytic energy, Fouvry 21 G9
    xp = rp/bc
    xa = ra/bc
    L0 = isochrone_L0(bc,M,astronomicalG)
    sp,sa = IsochroneSpSaFromRpRa(rp,ra,bc)
    return sqrt(2)*L0*xp*xa/sqrt((1+sp)*(1+sa)*(sp+sa))
end

"""
energy and angular momentum from isochrone model, using rpra
"""
function isochrone_EL_from_rpra(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    xp          = rp/bc
    xa          = ra/bc
    L0          = isochrone_L0(bc,M,astronomicalG)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    sp,sa       = IsochroneSpSaFromRpRa(rp,ra,bc)
    return scaleEnergy/(sp+sa),sqrt(2)*L0*xp*xa/sqrt((1+sp)*(1+sa)*(sp+sa))
end

"""
energy and angular momentum from isochrone model, using ae
"""
function IsochroneELFromAE(a::Float64,ecc::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    rp,ra = a*(1-ecc),a*(1+ecc)
    xp          = rp/bc
    xa          = ra/bc
    L0          = isochrone_L0(bc,M,astronomicalG)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    sp,sa       = IsochroneSpSaFromRpRa(rp,ra,bc)
    return scaleEnergy/(sp+sa),sqrt(2)*L0*xp*xa/sqrt((1+sp)*(1+sa)*(sp+sa))
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
    Omega1,Omega2 = IsochroneOmega12FromRpRa(rp,ra,bc,M,astronomicalG)
    sp,sa = IsochroneSpSaFromRpRa(rp,ra)
    return (3/sqrt(2))*(Omega1/Omega0)*(xr/sqrt(4-u^(2)))*(sqrt((sr+sp)*(sr+sa)*(sp+sa))/sqrt((xr+xp)*(xr+xa)))
end


"""
beta as a function of alpha for circular orbits
"""
function isochrone_beta_c(alpha::Float64)
    return 1/(1+alpha^(2/3))
end

"""
Theta function for the isochrone model
"""
function isochrone_drduINVvrfromrpra(rp::Float64,ra::Float64,u::Float64,bc::Float64=1.,Omega0::Float64=1.)
    Sigma, Delta = (ra+rp)*0.5, (ra-rp)*0.5 # Used for the mapping from u
    r = Sigma + Delta*henon_f(u) # Current value of the radius
    xp, xa, xr = rp/bc, ra/bc, r/bc # Rescaled pericentre, apocentre, and radius
    sqxp, sqxa, sqxr = sqrt(1.0+xp^(2)), sqrt(1.0+xa^(2)), sqrt(1.0+xr^(2)) # Pre-computing the values of sqrt(1+xp^2), sqrt(1+xa^2), and sqrt(1+xr^(2))
    #####
    drduINVvr = (3.0/(sqrt(2.0)))/(Omega0)*xr*sqrt(((sqxr+sqxp)*(sqxr+sqxa)*(sqxp+sqxa))/((xr+xp)*(xr+xa)*(4.0-u^(2))))# Analytical expression of (dr/du)(1/vr), that is always well-posed
    return drduINVvr # Output
end


"""
analytic jacobian for the change of variables from E,L to alpha,beta
"""
function IsochroneJacELtoAlphaBeta(alpha::Float64,
                                   beta::Float64,
                                   bc::Float64=1.,
                                   M::Float64=1.,
                                   astronomicalG::Float64=1.)

    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    scaleAction = isochrone_L0(bc,M,astronomicalG)
    Omega0      = isochrone_Omega0(bc,M,astronomicalG)

    return abs((1.0/6.0)*scaleEnergy*scaleAction/(alpha^(1/3)*(beta*(1.0-beta))^(3/2)))#*(1.0/Omega0)) # Output of the ABSOLUTE VALUE of the Jacobian. ATTENTION, contains the rescaling factor 1/Omega0
end


"""analytic definition of beta_c"""
function analytic_beta_c(x::Float64)::Float64
    return 1/(1 + x^(2/3))
end
