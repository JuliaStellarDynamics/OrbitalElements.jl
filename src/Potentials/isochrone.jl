"""
The isochrone potential definitions

  isochrone.jl is unique in that we can define many of the quantities analytically, which makes for a useful testbed for empirical potentials.

  nearly all quantities for the isochrone potential can be computed analytically

"""

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
function Ω₀Isochrone(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    return sqrt(astronomicalG*M/bc^3)
end

"""
isochrone energy scale, from Fouvry 21 (appendix G)
Emin = -GM/(2bc)
"""
function isochroneE0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #return -sqrt(astronomicalG*M/bc)
    return -astronomicalG*M/bc
end

"""
isochrone action scale, from Fouvry 21 (appendix G)
"""
function isochroneL0(bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
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
    E = isochroneEfromrpra(rp,ra,bc,M,astronomicalG)
    L = isochroneLfromrpra(rp,ra,bc,M,astronomicalG)
    return (astronomicalG*M/sqrt(-2E)) - 0.5 * (L + sqrt(L*L + 4*astronomicalG*M*bc))
end


function IsochroneActionsFromAE(a::Float64,ecc::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    rp,ra = a*(1-ecc),a*(1+ecc)
    xp          = rp/bc
    xa          = ra/bc
    L0          = isochroneL0(bc,M,astronomicalG)
    scaleEnergy = isochroneE0(bc,M,astronomicalG)
    sp,sa       = IsochroneSpSaFromRpRa(rp,ra,bc)
    Jval = IsochroneJrRpRa(rp,ra,bc,M,astronomicalG)
    return Jval,sqrt(2)*L0*xp*xa/sqrt((1+sp)*(1+sa)*(sp+sa))
end



"""
compute the dimensionless function for Omega1
(Fouvry 21 eq. G5)
"""
function IsochroneαRpRa(rp::Float64,ra::Float64,bc::Float64=1.)
    sp,sa = IsochroneSpSaFromRpRa(rp,ra,bc)
    return (2/(sp+sa))^(3/2)
end

"""
compute the dimensionless function for Omega1
(Fouvry 21 eq. G5)
"""
function IsochroneαAE(a::Float64,ecc::Float64,bc::Float64=1.)
    sp,sa = IsochroneSpSaFromRpRa(a*(1-ecc),a*(1+ecc),bc)
    return (2/(sp+sa))^(3/2)
end

"""
compute the dimensionless function for Omega2 from (rp,ra)
(Fouvry 21 eq. G7)
"""
function IsochroneβRpRa(rp::Float64,ra::Float64,bc::Float64=1.)
    xp = rp/bc
    xa = ra/bc
    sp,sa = IsochroneSpSaFromRpRa(rp,ra,bc)
    return (1/2)*(1+(xp*xa)/((1+sp)*(1+sa)))
end

"""
compute the dimensionless function for Omega2 from (a,e)
(Fouvry 21 eq. G7)
"""
function IsochroneβAE(a::Float64,ecc::Float64,bc::Float64=1.)
    xp = (a*(1-ecc))/bc
    xa = (a*(1+ecc))/bc
    sp,sa = IsochroneSpSaFromRpRa(a*(1-ecc),a*(1+ecc),bc)
    return (1/2)*(1+(xp*xa)/((1+sp)*(1+sa)))
end


"""
analytic function to return isochrone frequencies from (rp,ra)
"""
function IsochroneOmega12FromRpRa(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    Omega0   = Ω₀Isochrone(bc,M,astronomicalG)
    omegaae = IsochroneαRpRa(rp,ra,bc)
    etaae   = IsochroneβRpRa(rp,ra,bc)
    return omegaae*Omega0,omegaae*etaae*Omega0

end

"""
analytic function to return isochrone frequencies from (a,e)
"""
function IsochroneOmega12FromAE(a::Float64,ecc::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    Omega0   = Ω₀Isochrone(bc,M,astronomicalG)
    omegaae = IsochroneαAE(a,ecc,bc)
    etaae   = IsochroneβAE(a,ecc,bc)
    return omegaae*Omega0,omegaae*etaae*Omega0
end

"""
inversion of EL -> α,β function
"""
function isochroneELfromαβ(α::Float64,β::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochroneE0(bc,M,astronomicalG)
    scaleAction = isochroneL0(bc,M,astronomicalG)
    E = 0.5*scaleEnergy*(α)^(2/3) # Value of the energy
    L = scaleAction*(2.0*β-1.0)/(sqrt(β*(1.0-β))) # Value of the angular momentum
    return E, L # Output
end

"""

@IMPROVE, uses a floor to avoid any sqrt problems with circular orbits
"""
function isochronerprafromEL(E::Float64,L::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochroneE0(bc,M,astronomicalG)
    scaleAction = isochroneL0(bc,M,astronomicalG)
    xc = ((0.5*scaleEnergy)/E) - 1.0                                 # Value of xc
    eccov = sqrt(max(0,1.0 - (L/scaleAction)^(2)*(1.0/xc)*(1.0+(1.0/xc)))) # Value of overline{e}
    xp = sqrt(max(0,(2.0 + xc*(1.0-eccov))*(xc*(1.0-eccov))))             # Value of xp
    xa = sqrt(max(0,(2.0 + xc*(1.0+eccov))*(xc*(1.0+eccov))))             # Value of xa
    rp, ra = xp*bc, xa*bc                                            # Value of rp,ra
    return rp, ra # Output
end

"""
function to wrap (α,β)->(E,L)->(rp,ra)->(a,e) conversions for isochrone
"""
function IsochroneAEFromOmega1Omega2(omega1::Float64,omega2::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    Omega0= Ω₀Isochrone(bc,M,astronomicalG)
    E,L   = isochroneELfromαβ(omega1/Omega0,omega2/omega1,bc,M,astronomicalG)
    rp,ra = isochronerprafromEL(E,L,bc,M,astronomicalG)
    a,e   = AEfromRpRa(rp,ra)
    return a,e
end

"""
energy from isochrone model, using rpra
Fouvry 21 G9
"""
function isochroneEfromrpra(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochroneE0(bc,M,astronomicalG)
    sp,sa       = IsochroneSpSaFromRpRa(rp,ra,bc)
    return scaleEnergy/(sp+sa)
end

"""
angular momentum from isochrone model, using rpra
"""
function isochroneLfromrpra(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # isochrone analytic energy, Fouvry 21 G9
    xp = rp/bc
    xa = ra/bc
    L0 = isochroneL0(bc,M,astronomicalG)
    sp,sa = IsochroneSpSaFromRpRa(rp,ra,bc)
    return sqrt(2)*L0*xp*xa/sqrt((1+sp)*(1+sa)*(sp+sa))
end

"""
energy and angular momentum from isochrone model, using rpra
"""
function isochroneELfromrpra(rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    xp          = rp/bc
    xa          = ra/bc
    L0          = isochroneL0(bc,M,astronomicalG)
    scaleEnergy = isochroneE0(bc,M,astronomicalG)
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
    L0          = isochroneL0(bc,M,astronomicalG)
    scaleEnergy = isochroneE0(bc,M,astronomicalG)
    sp,sa       = IsochroneSpSaFromRpRa(rp,ra,bc)
    return scaleEnergy/(sp+sa),sqrt(2)*L0*xp*xa/sqrt((1+sp)*(1+sa)*(sp+sa))
end

function isochronedthetadufromrpra(r::Float64,u::Float64,rp::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #=

     the isochrone analytic Jacobian, Fouvry 21 G10
    =#
    xp = rp/bc
    xa = ra/bc
    xr = r/bc
    sr = sqrt(1+xr^(2))
    Omega0 = Ω₀Isochrone(bc,M,astronomicalG)
    Omega1,Omega2 = IsochroneOmega12FromRpRa(rp,ra,bc,M,astronomicalG)
    sp,sa = IsochroneSpSaFromRpRa(rp,ra)
    return (3/sqrt(2))*(Omega1/Omega0)*(xr/sqrt(4-u^(2)))*(sqrt((sr+sp)*(sr+sa)*(sp+sa))/sqrt((xr+xp)*(xr+xa)))
end

"""
Theta function for the isochrone model
"""
function isochronedrduINVvrfromrpra(rp::Float64,ra::Float64,u::Float64,bc::Float64=1.,Omega0::Float64=1.)
    a,e = AEFromRpRa(rp,ra)
    r = ru(u,a,e)
    xp, xa, xr = rp/bc, ra/bc, r/bc # Rescaled pericentre, apocentre, and radius
    sqxp, sqxa, sqxr = sqrt(1.0+xp^(2)), sqrt(1.0+xa^(2)), sqrt(1.0+xr^(2)) # Pre-computing the values of sqrt(1+xp^2), sqrt(1+xa^2), and sqrt(1+xr^(2))
    #####
    drduINVvr = (3.0/(sqrt(2.0)))/(Omega0)*xr*sqrt(((sqxr+sqxp)*(sqxr+sqxa)*(sqxp+sqxa))/((xr+xp)*(xr+xa)*(4.0-u^(2))))# Analytical expression of (dr/du)(1/vr), that is always well-posed
    return drduINVvr # Output
end


"""
analytic jacobian for the change of variables from E,L to α,β
"""
function IsochroneJacELtoαβ(α::Float64,
                                   β::Float64,
                                   bc::Float64=1.,
                                   M::Float64=1.,
                                   astronomicalG::Float64=1.)

    scaleEnergy = isochroneE0(bc,M,astronomicalG)
    scaleAction = isochroneL0(bc,M,astronomicalG)
    Omega0      = Ω₀Isochrone(bc,M,astronomicalG)

    return abs((1.0/6.0)*scaleEnergy*scaleAction/(α^(1/3)*(β*(1.0-β))^(3/2)))#*(1.0/Omega0)) # Output of the ABSOLUTE VALUE of the Jacobian. ATTENTION, contains the rescaling factor 1/Omega0
end


"""analytic definition of βc"""
function IsochroneβC(x::Float64)::Float64
    return 1/(1 + x^(2/3))
end

"""
the analytic expression for Theta for the isochrone profile
Fouvry, Hamilton, Rozier, Pichon eq. G10
this has the major advantage of always being well-posed for -1<u<1
"""
function ThetaRpRaIsochrone(rp::Float64,ra::Float64,
                            u::Float64;
                            bc::Float64=1.0,Ω₀::Float64=1.0)::Float64

    # compute helper quantities: used for the mapping from u
    a,e = AEFromRpRa(rp,ra)
    r = ru(u,a,e)

    # rescaled pericentre, apocentre, and radius
    xp, xa, xr = rp/bc, ra/bc, r/bc

    # pre-compute the values of sqrt(1+xp^2), sqrt(1+xa^2), and sqrt(1+xr^(2))
    sqxp, sqxa, sqxr = sqrt(1.0+xp^(2)), sqrt(1.0+xa^(2)), sqrt(1.0+xr^(2))

    # analytical expression of (dr/du)(1/vr), that is always well-posed
    drduINVvr = (3.0/(sqrt(2.0)))/(Ω₀)*xr*sqrt(((sqxr+sqxp)*(sqxr+sqxa)*(sqxp+sqxa))/((xr+xp)*(xr+xa)*(4.0-u^(2))))

    # output
    return drduINVvr

end




"""
# Returning the resonance frequency,
# along the circular orbit
"""
function OmgCircIsochrone(n1::Int64,n2::Int64,α::Float64)
    return n1*α + n2*α*IsochroneβC(α) # Output of the circular frequency
end


"""
# Function that computes omgmin and omgmax
# for a given resonance vector (n1,n2)
# ATTENTION, specific to the isochrone case
αcut would allow for controlling the minimum frequency probed for certain resonances, if interested.
"""
function FindWminWmaxIsochrone(n1::Int64,n2::Int64)
    αcut = 0.0 # Minimum value allowed for α
    omg1 = n1*αcut + n2*αcut*0.5                 # Bottom left value of omega
    omg2 = n1*αcut + n2*αcut*IsochroneβC(αcut) # Top left value of omega
    omg3 = n1           + n2          *0.5                 # Bottom right value of omega
    #####
    omgmin, omgmax = extrema((omg1,omg2,omg3)) # Computing the min/max values from the edge
    #####
    if ((n1 == 0) || (n2 == 0)) # In these cases, min/max are reached in the edges
        return omgmin, omgmax # Output
    else # In other cases, there might be a maximum along the circular orbits
        Delta = n2*(n2 - 24*n1) # Value of the discriminant of P(q)
        if (Delta <= 0.0) # There are no maximum along circular orbits
            return omgmin, omgmax
        else # The discriminant is positive, there might be an extremum along circular orbits
            qcut = 1.0 + αcut^(2/3) # Minimum value of q
            sqDelta = sqrt(Delta) # Square root of the discriminant
            qm = (-n2 - sqDelta)/(6.0*n1) # First  root of P(q)
            qp = (-n2 + sqDelta)/(6.0*n1) # Second root of P(q)
            #####
            qmin, qmax = minmax(qm,qp) # We sort the roots. At this stage, we necessarily have qmin <= 0.0
            #####
            if (qcut <= qmax <= 2.0) # We have reached an extremum along circular orbits in the allowed domain
                αcirc = (qmax - 1.0)^(3/2) # Value of α for which an extremum is reached along circular orbits
                omgcirc = n1*αcirc + n2*αcirc*IsochroneβC(αcirc) # Value of the extremum
                #####
                omgmin, omgmax = extrema((omgmin,omgmax,omgcirc)) # Determining the new minimum and maximum values of the resonance frequency
                return omgmin, omgmax # Output
            else # We do not reach an extremum along circular orbits in the allowed domain
                return omgmin, omgmax # Output
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
    αcut = 0.0
    omgmin, omgmax = FindWminWmaxIsochrone(n1,n2) # Getting omgmin, omgmax for the current resonance
    #####
    exph = 0.5*(omgmax + omgmin + u*(omgmax - omgmin)) # Expression of h(u)
    #####
    # For n2=0, the boundary takes a simple form
    if (n2 == 0)
        vm = 0.5 # Minimum bound
        vp = IsochroneβC(exph/n1) # Maximum bound
        return vm, vp # Output
    end
    #####
    # For n2!=0, the boundary is more complicated
    #####
    # First, we account for the easy constraints
    vm = αcut # Minimum bound
    vp = 1.0 # Maximum bound
    #####
    kappa = n1 + n2*0.5 # Definition of the constant kappa
    ####
    # Accounting for the third constraint
    if     (n2 < 0)
        if     (kappa < 0)
            vp = min(vp, exph/kappa) # Updating vp
        elseif (kappa > 0)
            vm = max(vm, exph/kappa) # Updating vm
        end
    #####
    elseif (n2 > 0)
    #####
        if     (kappa < 0)
            vm = max(vm, exph/kappa) # Updating vm
        elseif (kappa > 0)
            vp = min(vp, exph/kappa) # Updating vp
        end
    end
    #####
    # Accounting for the fourth constraint
    #####
    # The case n1=0 is easy to account for
    if (n1 == 0)
        ####
        omgl, omgr = (OmgCircIsochrone(n1,n2,αcut)-exph), (OmgCircIsochrone(n1,n2,1.0)-exph) # Values on the left/right of the interval
        ####
        if (omgl*omgr < 0.0) # We can search for a root
            vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exph,αcut,1.0) # Finding the transition coordinate
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
        omgl, omgr = (OmgCircIsochrone(n1,n2,αcut)-exph), (OmgCircIsochrone(n1,n2,1.0)-exph) # Values on the left/right of the interval
        ####
        if (omgl*omgr < 0.0) # We can search for a root
            vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exph,αcut,1.0) # Finding the transition coordinate
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
        qc = 1.0 + (αcut)^(2/3) # Left value of the range in q
        #####
        qm = (-n2 - sqrt(Delta))/(6.0*n1) # Value of q{-}
        qp = (-n2 + sqrt(Delta))/(6.0*n1) # Value of q{+}
        #####
        qmin, qmax = minmax(qm,qp) # Sorting the two roots
        #####
        # The considered range q \in [qc,2] is on the right of the two roots of P(q)
        if (qmax <= qc)
            #####
            omgl, omgr = (OmgCircIsochrone(n1,n2,αcut)-exph), (OmgCircIsochrone(n1,n2,1.0)-exph) # Values on the left/right of the interval
            #####
            if (omgl*omgr < 0.0) # We can search for a root
                vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exph,αcut,1.0) # Finding the transition coordinate
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
            omgl, omgr = (OmgCircIsochrone(n1,n2,αcut)-exph), (OmgCircIsochrone(n1,n2,1.0)-exph) # Values on the left/right of the interval
            #####
            if (omgl*omgr < 0.0) # We can search for a root
                vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exph,αcut,1.0) # Finding the transition coordinate
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
            omgLEFT = OmgCircIsochrone(n1,n2,αcut) # Value of OmgCircIsochrone at the left of the interaval
            omgMIDD = OmgCircIsochrone(n1,n2,vmax)      # Value of OmgCircIsochrone at the centre of the interval
            omgRGHT = OmgCircIsochrone(n1,n2,1.0)       # Value of OmgCircIsochrone at the right of the interval
            #####
            if     ((n1 < 0) && (n2 < 0))
                #####
                error("BUG -- getvmvp: CASE 1.")
                #####
            elseif ((n1 < 0) && (n2 > 0))
                #####
                if (omgLEFT <= exph <= omgMIDD)
                    vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exph,αcut,vmax) # Finding the transition
                    #####
                    vm = max(vm,vlim) # Updating vm
                end
                #####
                if (omgRGHT <= exph <= omgMIDD)
                    vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exph,vmax,1.0) # Finding the transition coordinate
                    #####
                    vp = min(vp,vlim) # Updating vp
                end
                #####
            elseif ((n1 > 0) && (n2 < 0))
                #####
                if (omgMIDD <= exph <= omgLEFT)
                    vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exph,αcut,vmax) # Finding the transition coordinate
                    #####
                    vm = max(vm,vlim) # Updating vm
                end
                #####
                if (omgMIDD <= exph <= omgRGHT)
                    vlim = bisection(v -> OmgCircIsochrone(n1,n2,v)-exph,vmax,1.0) # Finding the transition coordinate
                    #####
                    vp = min(vp,vlim) # Updating vp
                end
                #####
            elseif ((n1 > 0) && (n2 > 0))
                #####
                error("BUG -- getvmvp: Case 2.")
            end
            #####
            return vm, vp # Output
        end
    end
end


function GetϖIsochrone(omg::Complex{Float64},
                           n1::Int64,n2::Int64)

    wmin,wmax = FindWminWmaxIsochrone(n1,n2)

    return Getϖ(omg,wmin,wmax)

end
