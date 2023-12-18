"""
The isochrone potential definitions

  isochrone.jl is unique in that we can define many of the quantities analytically, which makes for a useful testbed for empirical potentials.

  nearly all quantities for the isochrone potential can be computed analytically

"""

#####################################
# Isochrone structures
#####################################
"""
Abstract Isochrone potential structure
"""
abstract type IsochronePotential <: CentralCorePotential end

"""
Isochrone potential structure using numerical computations
"""
struct NumericalIsochrone <: IsochronePotential
    G::Float64      # Gravitational constant
    M::Float64      # Total mass
    bc::Float64     # Characteristic radius
end

"""
Isochrone potential structure using analytical computations
"""
struct AnalyticalIsochrone <: IsochronePotential
    G::Float64      # Gravitational constant
    M::Float64      # Total mass
    bc::Float64     # Characteristic radius
end

"""
    NumericalIsochrone([, bc, M, G])

Create an Isochrone potential structure with characteristic radius `bc`
total mass `M` and gravitational constant `G`.
"""
function NumericalIsochrone(;G::Float64=1.,M::Float64=1.,bc::Float64=1.)
    return NumericalIsochrone(G,M,bc)
end

"""
    AnalyticalIsochrone([, bc, M, G])

Create an Isochrone potential structure with characteristic radius `bc`
total mass `M` and gravitational constant `G`.

For this structure, most computations will use analytical expressions.
"""
function AnalyticalIsochrone(;G::Float64=1.,M::Float64=1.,bc::Float64=1.)
    return AnalyticalIsochrone(G,M,bc)
end

#####################################
# Potential methods for the Isochrone
#####################################
function ψ(model::IsochronePotential,r::Float64)

    x = r/model.bc
    scale = model.G*model.M/model.bc

    return - scale / (1.0+sqrt(1.0+x^2))
end

function dψ(model::IsochronePotential,r::Float64)
    
    x = r/model.bc
    scale = model.G*model.M/(model.bc^2)

    return scale / (sqrt(1.0+x^(-2))*(1.0+sqrt(1.0+x^2))^2)
end

function d2ψ(model::IsochronePotential,r::Float64)

    x = r/model.bc
    scale = model.G*model.M/(model.bc^3)
    sqxp = 1.0 + x^2
    invx = x^(-1)
    invsqxp = 1.0 + invx^2
    return scale *(- ((sqxp^(3/2))*(invx+sqrt(invsqxp))^2)^(-1)
                    - 2*(invsqxp*(1.0+sqrt(sqxp))^3)^(-1)
                    + (sqrt(sqxp)*(1.0+sqrt(sqxp))^2)^(-1))
end

#####################################
# Scales for the Isochrone
#####################################
function Ω₀(model::IsochronePotential)
    return sqrt(model.G*model.M/(model.bc^3))
end

"""
    E₀(model::IsochronePotential)

for isochrone, see Fouvry 21 (appendix G).
"""
function E₀(model::IsochronePotential)
    return -model.G*model.M/model.bc
end

"""
    L₀(model::IsochronePotential)

for isochrone, see Fouvry 21 (appendix G)
"""
function L₀(model::IsochronePotential)
    return sqrt(model.G*model.M*model.bc)
end


#####################################
# Analytical functions for Isochrone
#####################################

"""
    SpSaFromRpRa(model::AnalyticalIsochrone, rp, ra)

isochrone reduced coordinates from pericentre `rp` and apocentre `ra`.
"""
function SpSaFromRpRa(model::AnalyticalIsochrone,rp::Float64,ra::Float64)
    xp, xa = (rp, ra) ./ model.bc
    return sqrt(1+xp^2), sqrt(1+xa^2)
end

"""
    SpSafromAE(model::AnalyticalIsochrone, a, e)

isochrone reduced coordinates from semi-major axis `a` and eccentricity `e`.
"""
function SpSaFromAE(model::AnalyticalIsochrone,a::Float64,e::Float64)
    rp, ra = RpRaFromAE(a,e)
    return SpSaFromRpRa(model,rp,ra)
end

"""   
    EFromAE(model::AnalyticalIsochrone, a, e)

for isochrone analytical version, see equation (G9)
in Fouvry&Prunet (2021)
"""
function EFromAE(model::AnalyticalIsochrone,a::Float64,e::Float64)
    sp,sa = SpSaFromAE(model,a,e)
    return E₀(model)/(sp+sa)
end

"""   
    LFromAE(model::AnalyticalIsochrone, a, e)

for isochrone analytical version, see equation (G9)
in Fouvry&Prunet (2021)
"""
function LFromAE(model::AnalyticalIsochrone,a::Float64,e::Float64)

    rp, ra = RpRaFromAE(a,e)
    xp, xa = (rp, ra) ./ model.bc
    sp,sa = SpSaFromRpRa(model,rp,ra)

    return sqrt(2)*L₀(model)*xp*xa/sqrt((1+sp)*(1+sa)*(sp+sa))
end

"""
    ELFromAE(model::AnalyticalIsochrone, a, e)

for isochrone analytical version, see equation (G9)
in Fouvry&Prunet (2021)
"""
function ELfromAE(model::AnalyticalIsochrone,a::Float64,e::Float64)
    return EFromAE(model,a,e), LFromAE(model,a,e)
end

"""
    JFromAE(model::AnalyticalIsochrone, a, e)

for isochrone analytical version, see equation (G3)
in Fouvry&Prunet (2021)
"""
function JFromAE(model::AnalyticalIsochrone,a::Float64,e::Float64)

    E, L = ELFromAE(model,a,e)
    return (model.G*model.M/sqrt(-2E)) - 0.5 * (L + sqrt(L*L + 4*model.G*model.M*model.bc))
end

"""
    ComputeActionsAE(model::AnalyticalIsochrone, a, e)

for isochrone analytical version, see equation (G3)
in Fouvry&Prunet (2021)
"""
function ComputeActionsAE(model::AnalyticalIsochrone,a::Float64,e::Float64)

    Jval = JFromAE(model,a,e)
    return Jval, LFromAE(model,a,e)
end

"""
    αβFromAE(model::AnalyticalIsochrone, a, e)

for isochrone analytical version, 
see equations (G5-G7) in Fouvry&Prunet (2021)
"""
function αβFromAE(model::AnalyticalIsochrone,a::Float64,e::Float64)

    rp, ra = RpRaFromAE(a,e)
    xp, xa = (rp, ra) ./ model.bc
    sp,sa = SpSaFromRpRa(model,rp,ra)

    return (2/(sp+sa))^(3/2), (1/2)*(1+(xp*xa)/((1+sp)*(1+sa)))
end

function ComputeFrequenciesAE(model::AnalyticalIsochrone,
                              a::Float64,e::Float64)
    
    α, β = αβFromAE(model,a,e)
    return FrequenciesFromαβ(α,β,Ω₀(model))
end

"""
    ELFromαβ(model::AnalyticalIsochrone, α, β)

for isochrone analytical version, 
see equations @TOCOMPLETE in Fouvry&Prunet (2021)
"""
function ELFromαβ(model::AnalyticalIsochrone,α::Float64,β::Float64)

    return 0.5*E₀(model)*(α)^(2/3), L₀(model)*(2.0*β-1.0)/(sqrt(β*(1.0-β)))
end

"""    
    AEFromEL(model::AnalyticalIsochrone, E, L)

for isochrone analytical version, 
see equations @TOCOMPLETE in Fouvry&Prunet (2021)

@IMPROVE, uses a floor to avoid any sqrt problems with circular orbits
"""
function AEFromEL(model::AnalyticalIsochrone,E::Float64,L::Float64)

    xc = (0.5*E₀(model)/E) - 1.0
    eccov = sqrt(max(0,1.0 - (L/L₀(model))^(2)*(1.0/xc)*(1.0+(1.0/xc))))
    xp = sqrt(max(0.,(2.0 + xc*(1.0-eccov))*(xc*(1.0-eccov))))
    xa = sqrt(max(0.,(2.0 + xc*(1.0+eccov))*(xc*(1.0+eccov))))
    rp, ra = xp*model.bc, xa*model.bc

    return AEFromRpRa(rp,ra)
end

function ComputeAEFromFrequencies(model::AnalyticalIsochrone,Ω1::Float64,Ω2::Float64)
    
    # (Ω1,Ω2) ↦ (α,β), generically analytical
    α, β = αβFromFrequencies(Ω1,Ω2,Ω₀(model))
    # (α,β) ↦ (E,L), analytic expressions specific to isochrone potential
    E, L = ELFromαβ(model,α,β)
    # (E,L) ↦ (a,e), analytic expressions specific to isochrone potential
    a, e = AEFromEL(model,E,L)

    return a,e
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
