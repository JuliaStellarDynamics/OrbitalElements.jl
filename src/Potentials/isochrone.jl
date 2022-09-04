#=
The isochrone potential definitions

  isochrone.jl is unique in that we can define many of the quantities analytically, which makes for a useful testbed for empirical potentials.

  nearly all quantities for the isochrone potential can be computed analytically

=#

"""ψIsochrone(r[,bc,M,G])

the isochrone potential
"""
function ψIsochrone(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    rbc = r^2 + bc^2
    return -astronomicalG*M*(bc+sqrt(rbc))^(-1)
end

"""dψIsochrone(r[,bc,M,G])

the isochrone potential derivative
"""
function dψIsochrone(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    rbc = r^2 + bc^2
    return astronomicalG*M*r*(sqrt(rbc)*(sqrt(rbc)+bc)^2)^(-1)
end

"""d2ψIsochrone(r[,bc,M,G])

the isochrone potential second derivative
"""
function d2ψIsochrone(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=
    =#
    rbc = r^2 + bc^2
    return astronomicalG*M*(- r^(2)*((rbc^(3/2))*(sqrt(rbc)+bc)^2)^(-1)
                          - 2*r^(2)*(rbc*(sqrt(rbc)+bc)^3)^(-1)
                          + (sqrt(rbc)*(sqrt(rbc)+bc)^2)^(-1))
end

"""d3ψIsochrone(r[,bc,M,G])

the isochrone potential third derivative
"""
function d3ψIsochrone(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
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

"""d4ψIsochrone(r[,bc,M,G])

the isochrone potential fourth derivative
"""
function d4ψIsochrone(r::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
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
function Omega0Isochrone(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
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
function IsochroneJrRpRa(rp::Float64,ra::Float64,
                         bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
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
    Omega0   = Omega0Isochrone(bc,M,astronomicalG)
    omega_ae = IsochroneAlphaRpRa(rp,ra,bc)
    eta_ae   = IsochroneBetaRpRa(rp,ra,bc)
    return omega_ae*Omega0,omega_ae*eta_ae*Omega0

end

"""
analytic function to return isochrone frequencies from (a,e)
"""
function IsochroneOmega12FromAE(a::Float64,ecc::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    Omega0   = Omega0Isochrone(bc,M,astronomicalG)
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
function IsochroneAEFromOmega1Omega2(omega1::Float64,omega2::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    Omega0= Omega0Isochrone(bc,M,astronomicalG)
    E,L   = isochrone_EL_from_alphabeta(omega1/Omega0,omega2/omega1,bc,M,astronomicalG)
    rp,ra = isochrone_rpra_fromEL(E,L,bc,M,astronomicalG)
    a,e   = AEfromRpRa(rp,ra)
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
    Omega0 = Omega0Isochrone(bc,M,astronomicalG)
    Omega1,Omega2 = IsochroneOmega12FromRpRa(rp,ra,bc,M,astronomicalG)
    sp,sa = IsochroneSpSaFromRpRa(rp,ra)
    return (3/sqrt(2))*(Omega1/Omega0)*(xr/sqrt(4-u^(2)))*(sqrt((sr+sp)*(sr+sa)*(sp+sa))/sqrt((xr+xp)*(xr+xa)))
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
    Omega0      = Omega0Isochrone(bc,M,astronomicalG)

    return abs((1.0/6.0)*scaleEnergy*scaleAction/(alpha^(1/3)*(beta*(1.0-beta))^(3/2)))#*(1.0/Omega0)) # Output of the ABSOLUTE VALUE of the Jacobian. ATTENTION, contains the rescaling factor 1/Omega0
end


"""analytic definition of beta_c"""
function IsochroneBetaC(x::Float64)::Float64
    return 1/(1 + x^(2/3))
end

"""
the analytic expression for Theta for the isochrone profile
Fouvry, Hamilton, Rozier, Pichon eq. G10
this has the major advantage of always being well-posed for -1<u<1
"""
function ThetaRpRaIsochrone(rp::Float64,ra::Float64,
                            u::Float64;
                            bc::Float64=1.0,Ω0::Float64=1.0)::Float64

    # compute helper quantities: used for the mapping from u
    Sigma, Delta = (ra+rp)*0.5, (ra-rp)*0.5

    # Current value of the radius
    r = Sigma + Delta*henon_f(u)

    # rescaled pericentre, apocentre, and radius
    xp, xa, xr = rp/bc, ra/bc, r/bc

    # pre-compute the values of sqrt(1+xp^2), sqrt(1+xa^2), and sqrt(1+xr^(2))
    sqxp, sqxa, sqxr = sqrt(1.0+xp^(2)), sqrt(1.0+xa^(2)), sqrt(1.0+xr^(2))

    # analytical expression of (dr/du)(1/vr), that is always well-posed
    drduINVvr = (3.0/(sqrt(2.0)))/(Ω0)*xr*sqrt(((sqxr+sqxp)*(sqxr+sqxa)*(sqxp+sqxa))/((xr+xp)*(xr+xa)*(4.0-u^(2))))

    # output
    return drduINVvr

end




"""
# Returning the resonance frequency,
# along the circular orbit
"""
function OmgCircIsochrone(n1::Int64,n2::Int64,alpha::Float64)
    return n1*alpha + n2*alpha*IsochroneBetaC(alpha) # Output of the circular frequency
end


"""
# Function that computes omg_min and omg_max
# for a given resonance vector (n1,n2)
# ATTENTION, specific to the isochrone case
alpha_cut would allow for controlling the minimum frequency probed for certain resonances, if interested.
"""
function FindWminWmaxIsochrone(n1::Int64,n2::Int64)
    alpha_cut = 0.0 # Minimum value allowed for alpha
    omg_1 = n1*alpha_cut + n2*alpha_cut*0.5                 # Bottom left value of omega
    omg_2 = n1*alpha_cut + n2*alpha_cut*IsochroneBetaC(alpha_cut) # Top left value of omega
    omg_3 = n1           + n2          *0.5                 # Bottom right value of omega
    #####
    omg_min, omg_max = extrema((omg_1,omg_2,omg_3)) # Computing the min/max values from the edge
    #####
    if ((n1 == 0) || (n2 == 0)) # In these cases, min/max are reached in the edges
        return omg_min, omg_max # Output
    else # In other cases, there might be a maximum along the circular orbits
        Delta = n2*(n2 - 24*n1) # Value of the discriminant of P(q)
        if (Delta <= 0.0) # There are no maximum along circular orbits
            return omg_min, omg_max
        else # The discriminant is positive, there might be an extremum along circular orbits
            q_cut = 1.0 + alpha_cut^(2/3) # Minimum value of q
            sqDelta = sqrt(Delta) # Square root of the discriminant
            qm = (-n2 - sqDelta)/(6.0*n1) # First  root of P(q)
            qp = (-n2 + sqDelta)/(6.0*n1) # Second root of P(q)
            #####
            qmin, qmax = minmax(qm,qp) # We sort the roots. At this stage, we necessarily have qmin <= 0.0
            #####
            if (q_cut <= qmax <= 2.0) # We have reached an extremum along circular orbits in the allowed domain
                alpha_circ = (qmax - 1.0)^(3/2) # Value of alpha for which an extremum is reached along circular orbits
                omg_circ = n1*alpha_circ + n2*alpha_circ*IsochroneBetaC(alpha_circ) # Value of the extremum
                #####
                omg_min, omg_max = extrema((omg_min,omg_max,omg_circ)) # Determining the new minimum and maximum values of the resonance frequency
                return omg_min, omg_max # Output
            else # We do not reach an extremum along circular orbits in the allowed domain
                return omg_min, omg_max # Output
            end
        end
    end
end


"""
# Function that returns the edges (vm,vp)
# of the integration over v,
# for a given resonance (n1,n2) and coordinate u
"""
function FindVminVmaxIsochrone(n1::Int64,n2::Int64,u::Float64)
    alpha_cut = 0.0
    omg_min, omg_max = FindWminWmaxIsochrone(n1,n2) # Getting omg_min, omg_max for the current resonance
    #####
    exp_h = 0.5*(omg_max + omg_min + u*(omg_max - omg_min)) # Expression of h(u)
    #####
    # For n2=0, the boundary takes a simple form
    if (n2 == 0)
        vm = 0.5 # Minimum bound
        vp = IsochroneBetaC(exp_h/n1) # Maximum bound
        return vm, vp # Output
    end
    #####
    # For n2!=0, the boundary is more complicated
    #####
    # First, we account for the easy constraints
    vm = alpha_cut # Minimum bound
    vp = 1.0 # Maximum bound
    #####
    kappa = n1 + n2*0.5 # Definition of the constant kappa
    ####
    # Accounting for the third constraint
    if     (n2 < 0)
        if     (kappa < 0)
            vp = min(vp, exp_h/kappa) # Updating vp
        elseif (kappa > 0)
            vm = max(vm, exp_h/kappa) # Updating vm
        end
    #####
    elseif (n2 > 0)
    #####
        if     (kappa < 0)
            vm = max(vm, exp_h/kappa) # Updating vm
        elseif (kappa > 0)
            vp = min(vp, exp_h/kappa) # Updating vp
        end
    end
    #####
    # Accounting for the fourth constraint
    #####
    # The case n1=0 is easy to account for
    if (n1 == 0)
        ####
        omg_l, omg_r = (OmgCircIsochrone(n1,n2,alpha_cut)-exp_h), (OmgCircIsochrone(n1,n2,1.0)-exp_h) # Values on the left/right of the interval
        ####
        if (omg_l*omg_r < 0.0) # We can search for a root
            vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exp_h,alpha_cut,1.0) # Finding the transition coordinate
            vm = max(vm,vlim) # Updating vm
        end
        #####
        return vm, vp # Output
    end
    #####
    # Now, we know for sure that n1!=0 and n2!=0
    Delta = n2*(n2-24*n1) # Value of the discriminant
    #####
    if (Delta <= 0) # There is no extrema along circular orbits
        ####
        omg_l, omg_r = (OmgCircIsochrone(n1,n2,alpha_cut)-exp_h), (OmgCircIsochrone(n1,n2,1.0)-exp_h) # Values on the left/right of the interval
        ####
        if (omg_l*omg_r < 0.0) # We can search for a root
            vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exp_h,alpha_cut,1.0) # Finding the transition coordinate
            #####
            if     (n1*n2 > 0)
                vm = max(vm,vlim) # Updating vm
            elseif (n1*n2 < 0)
                vp = min(vp,vlim) # Updating vp
            end
            ######
        end
        ####
        return vm, vp # Output
    #####
    # Otherwise, we have Delta > 0, and we might find an extremum along circular orbits
    else
        qc = 1.0 + (alpha_cut)^(2/3) # Left value of the range in q
        #####
        qm = (-n2 - sqrt(Delta))/(6.0*n1) # Value of q_{-}
        qp = (-n2 + sqrt(Delta))/(6.0*n1) # Value of q_{+}
        #####
        qmin, qmax = minmax(qm,qp) # Sorting the two roots
        #####
        # The considered range q \in [qc,2] is on the right of the two roots of P(q)
        if (qmax <= qc)
            #####
            omg_l, omg_r = (OmgCircIsochrone(n1,n2,alpha_cut)-exp_h), (OmgCircIsochrone(n1,n2,1.0)-exp_h) # Values on the left/right of the interval
            #####
            if (omg_l*omg_r < 0.0) # We can search for a root
                vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exp_h,alpha_cut,1.0) # Finding the transition coordinate
                #####
                if     (n1*n2 > 0)
                    vm = max(vm,vlim) # Updating vm
                elseif (n1*n2 < 0)
                    vp = min(vp,vlim) # Updating vp
                end
                #####
            end
            #####
            return vm, vp # Output
        #####
        # The considered range q \in [qc,2] is within the two roots of P(q)
        # The only change with the above case is a change in the sign
        # for the condition on n1*n2
        elseif (2.0 <= qmax)
            #####
            omg_l, omg_r = (OmgCircIsochrone(n1,n2,alpha_cut)-exp_h), (OmgCircIsochrone(n1,n2,1.0)-exp_h) # Values on the left/right of the interval
            #####
            if (omg_l*omg_r < 0.0) # We can search for a root
                vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exp_h,alpha_cut,1.0) # Finding the transition coordinate
                #####
                if     (n1*n2 > 0)
                    vp = min(vp,vlim) # Updating vp
                elseif (n1*n2 < 0)
                    vm = max(vm,vlim) # Updating vm
                end
                #####
            end
            ####
            return vm, vp # Output
        #####
        # The considered range q \in [qc,2] contains qmax,
        # i.e. contains a root of P(q)
        # This is the case that is more difficult to consider
        else
            vmax = (qmax - 1.0)^(3/2) # Position in v where the monotonicity changes
            #####
            omg_LEFT = OmgCircIsochrone(n1,n2,alpha_cut) # Value of OmgCircIsochrone at the left of the interaval
            omg_MIDD = OmgCircIsochrone(n1,n2,vmax)      # Value of OmgCircIsochrone at the centre of the interval
            omg_RGHT = OmgCircIsochrone(n1,n2,1.0)       # Value of OmgCircIsochrone at the right of the interval
            #####
            if     ((n1 < 0) && (n2 < 0))
                #####
                error("BUG -- get_vmvp: CASE 1.")
                #####
            elseif ((n1 < 0) && (n2 > 0))
                #####
                if (omg_LEFT <= exp_h <= omg_MIDD)
                    vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exp_h,alpha_cut,vmax) # Finding the transition
                    #####
                    vm = max(vm,vlim) # Updating vm
                end
                #####
                if (omg_RGHT <= exp_h <= omg_MIDD)
                    vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exp_h,vmax,1.0) # Finding the transition coordinate
                    #####
                    vp = min(vp,vlim) # Updating vp
                end
                #####
            elseif ((n1 > 0) && (n2 < 0))
                #####
                if (omg_MIDD <= exp_h <= omg_LEFT)
                    vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exp_h,alpha_cut,vmax) # Finding the transition coordinate
                    #####
                    vm = max(vm,vlim) # Updating vm
                end
                #####
                if (omg_MIDD <= exp_h <= omg_RGHT)
                    vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exp_h,vmax,1.0) # Finding the transition coordinate
                    #####
                    vp = min(vp,vlim) # Updating vp
                end
                #####
            elseif ((n1 > 0) && (n2 > 0))
                #####
                error("BUG -- get_vmvp: Case 2.")
            end
            #####
            return vm, vp # Output
        end
    end
end


function GetVarpiIsochrone(omg::Complex{Float64},
                           n1::Int64,n2::Int64)

    w_min,w_max = FindWminWmaxIsochrone(n1,n2)

    return GetVarpi(omg,w_min,w_max)

end
