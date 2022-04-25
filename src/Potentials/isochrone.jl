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
function isochrone_spsa_from_rpra(rp::Float64,ra::Float64,bc::Float64=1.)
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
function isochrone_omega_ae(rp::Float64,ra::Float64,bc::Float64=1.)
    sp,sa = isochrone_spsa_from_rpra(rp,ra,bc)
    return (2/(sp+sa))^(3/2)
end

"""
compute the dimensionless function for Omega2
(Fouvry 21 eq. G7)
"""
function isochrone_eta_ae(rp::Float64,ra::Float64,bc::Float64=1.)
    xp = rp/bc
    xa = ra/bc
    sp,sa = isochrone_spsa_from_rpra(rp,ra,bc)
    return (1/2)*(1+(xp*xa)/((1+sp)*(1+sa)))
end

"""
analytic function to return isochrone frequencies
"""
function isochrone_Omega_1_2(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    Omega0   = isochrone_Omega0(bc,M,astronomicalG)
    omega_ae = isochrone_omega_ae(rp,ra,bc)
    eta_ae   = isochrone_eta_ae(rp,ra,bc)
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
    sp,sa       = isochrone_spsa_from_rpra(rp,ra,bc)
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
    sp,sa = isochrone_spsa_from_rpra(rp,ra,bc)
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
    sp,sa       = isochrone_spsa_from_rpra(rp,ra,bc)
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
    Omega1,Omega2 = isochrone_Omega_1_2(rp,ra,bc,M,astronomicalG)
    sp,sa = isochrone_spsa_from_rpra(rp,ra)
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
this is for the change of variables from E,L to alpha,beta
"""
function isochrone_JacEL_to_alphabeta(alpha::Float64,beta::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    scaleAction = isochrone_L0(bc,M,astronomicalG)
    Omega0      = isochrone_Omega0(bc,M,astronomicalG)
    return abs((1.0/6.0)*scaleEnergy*scaleAction/(alpha^(1/3)*(beta*(1.0-beta))^(3/2)))#*(1.0/Omega0)) # Output of the ABSOLUTE VALUE of the Jacobian. ATTENTION, contains the rescaling factor 1/Omega0
end


"""
"""
function isochrone_Q_ROI(E::Float64,L::Float64,Ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    Q = (E + L^(2)/(2.0*Ra^(2)))/scaleEnergy # Computing the Q parameter
    return Q # Output
end


"""
Saha distribution function
ra is the anisotropy radius
"""
function isochrone_Saha_DF(E::Float64,L::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #scaleEnergy = isochrone_E0(bc,M,astronomicalG)

    Q = isochrone_Q_ROI(E,L,ra,bc,M,astronomicalG)
    gamma = (bc/ra)^2

    prefactor = (M/((G*M*bc)^(3/2))) * (1/(128*sqrt(2)*pi)) * (sqrt(Q)/((1-Q)^4))

    term1 = 27 + 77gamma - (66+286gamma)*Q + (320+136gamma)*Q*Q -(240+32gamma)*Q*Q*Q + 64*Q*Q*Q*Q
    term2 = ((3*asin(sqrt(Q)))/sqrt(Q*(1-Q))) * ((-9 + 17gamma) + (28-44gamma)*Q + (16-8gamma)*Q*Q)

    return prefactor * (term1+term2)
end


"""
we also need dF/dQ and then we have the isochrone complete!
@ATTENTION, check the pi^3 factor here
@ATTENTION, check this answer against isochrone_Saha_DF
"""
function DFISOQ_ROI(Q::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    gamma_ROI = (bc/ra)^2
    return (Mtot*sqrt(Q)*(27.0+77.0*gamma_ROI+2.0*Q*(-11.0*(3.0+13.0*gamma_ROI)+
    4.0*Q*(40+2.0*Q*(-15.0+4.0*Q-2.0*gamma_ROI)+17.0*gamma_ROI))+
    (3.0*(-9.0+4.0*Q*(7.0-2.0*Q*(-2.0+gamma_ROI)-11.0*gamma_ROI)+17.0*gamma_ROI)*
    asin(sqrt(Q)))/sqrt(-((-1.0+Q)*Q))))/(128.0*sqrt(2.0)*(bc*astronomicalG*M)^(3/2)*(pi)^(3)*(-1.0+Q)^(4))
end

"""
Saha distribution function derivative w.r.t. Q
ra is the anisotropy radius
"""
function isochrone_Saha_dDFdQ(Q::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    gamma_ROI = (bc/ra)^2
    return (M*((-1.0+Q)*(-128.0*gamma_ROI+Q*(-75.0+451.0*gamma_ROI+2.0*Q*
    (-659.0+387.0*gamma_ROI+4.0*Q*(90-23.0*gamma_ROI+2.0*Q*(-21.0+4.0*Q+2.0*gamma_ROI)))))-
    15.0*sqrt(-((-1.0+Q)*Q))*(5.0-13.0*gamma_ROI+4.0*Q*(-13.0+2.0*Q*(-2.0+gamma_ROI)+
    17.0*gamma_ROI))*asin(sqrt(Q))))/(256.0*sqrt(2.0)*(bc*astronomicalG*M)^(3/2)*(pi)^(3)*(-1.0+Q)^(6)*sqrt(Q))
end

"""
jacobian for converting dF/dQ dQ/dE -> dF/dE
"""
function isochrone_dQdE_ROI(E::Float64,L::Float64,Ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    return 1.0/scaleEnergy # Output
end


"""
jacobian for converting dF/dQ dQ/dL -> dF/dL
"""
function isochrone_dQdL_ROI(E::Float64,L::Float64,Ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    return (L)/(scaleEnergy*Ra^(2)) # Output
end




"""isochrone_isotropic_DF(E[,bc,M,G])
the isotropic DF
"""
function isochrone_isotropic_DF(E::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    mE = E/scaleEnergy # Dimensionless energy that goes between 0 and 1/2
    return (sqrt(mE)*M*(27.0+2.0*mE*(-1.0+4.0*mE)*(33.0+4.0*mE*(-7.0+2.0*mE))+
    (3.0*(-9.0+4.0*mE*(7.0+4.0*mE))*asin(sqrt(mE)))/sqrt(-((-1.0+mE)*mE))))/
    (128.0*sqrt(2.0)*(-1.0+mE)^(4)*(bc*astronomicalG*M)^(3/2)*(pi)^(3))
end

"""isochrone_isotropic_dDFdE(E[,bc,M,G])
the isotropic DF derivative w.r.t. E
"""
function isochrone_isotropic_dDFdE(E::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochrone_E0(bc,M,astronomicalG)
    mE = E/scaleEnergy # Dimensionless energy that goes between 0 and 1/2
    return (M*((-1.0+mE)*sqrt(mE)*(-75.0+2.0*mE*(-659.0+8.0*mE*(45.0+mE*(-21.0+4.0*mE))))+
    15.0*sqrt(1.0-mE)*(-5.0+4.0*mE*(13.0+4.0*mE))*asin(sqrt(mE))))/
    (256.0*sqrt(2.0)*scaleEnergy*(-1.0+mE)^(6)*(bc*astronomicalG*M)^(3/2)*(pi)^(3))
end
