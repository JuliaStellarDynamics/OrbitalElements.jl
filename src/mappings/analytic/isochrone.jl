
########################################################################
#
# Appropriate reduced coordinates
#
########################################################################
"""
    spsa_from_rpra(rp, ra, model::AnalyticIsochrone)

isochrone reduced coordinates from pericentre `rp` and apocentre `ra`.
"""
function spsa_from_rpra(rp::Float64, ra::Float64, model::AnalyticIsochrone)
    xp, xa = (rp, ra) ./ model.bc
    return sqrt(1 + xp^2), sqrt(1 + xa^2)
end

"""
    spsa_from_ae(a, e, model::AnalyticIsochrone)

isochrone reduced coordinates from semi-major axis `a` and eccentricity `e`.
"""
function spsa_from_ae(a::Float64, e::Float64, model::AnalyticIsochrone)
    rp, ra = rpra_from_ae(a, e)
    return spsa_from_rpra(rp,ra, model)
end

########################################################################
#
# (a,e) ↔ (E,L) mappings
#
########################################################################
"""
for isochrone analytical version, see equation (G9) in Fouvry & Prunet (2022)
"""
function EL_from_ae(
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    rp, ra = rpra_from_ae(a, e)
    xp, xa = (rp, ra) ./ model.bc
    sp, sa = spsa_from_rpra(rp, ra, model)
    E = energy_scale(model) / (sp + sa)
    L = sqrt(2) * momentum_scale(model) * xp * xa / sqrt((1 + sp) * (1 + sa) * (sp + sa))
    return E, L 
end


"""    
for isochrone analytical version, 
see equations @TOCOMPLETE in Fouvry&Prunet (2021)
@IMPROVE, uses a floor to avoid any sqrt problems with circular orbits
"""
function ae_from_EL(
    E::Float64,
    L::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    xc = energy_scale(model) / (2E) - 1
    eccov = sqrt(max(0.0, 1 - (L / momentum_scale(model))^2 * (1 / xc) * (1 + (1 / xc))))
    xp = sqrt(max(0.0,(2 + xc * (1 - eccov)) * (xc * (1 - eccov))))
    xa = sqrt(max(0.0,(2 + xc * (1 + eccov)) * (xc * (1 + eccov))))
    rp, ra = xp * model.bc, xa * model.bc
    return ae_from_rpra(rp, ra)
end

########################################################################
#
# (a,e) ↔ (J,L) mappings
#
########################################################################
"""
for isochrone analytical version, see equation (G3) in Fouvry & Prunet (2022)
"""
function _radial_action_from_ae(
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    E, L = EL_from_ae(a, e, model, params)
    G, M, bc = model.G, model.M, model.bc
    return G * M / sqrt(-2E) - (L + sqrt(L^2 + 4G * M * bc)) / 2
end

########################################################################
#
# (a,e) ↔ frequencies mappings: integrands
#
########################################################################
"""
the analytic expression for Theta for the isochrone profile are given in 
Fouvry, Hamilton, Rozier, Pichon eq. (G10). It has the major advantage of 
always being well-posed for -1<u<1
"""
function Θ(
    u::Float64,
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    bc = model.bc
    rp, ra = rpra_from_ae(a, e)
    r = radius_from_anomaly(u,a,e)

    xp, xa, xr = rp/bc, ra/bc, r/bc # Rescaled pericentre, apocentre, and radius
    # Pre-computing
    sqxp, sqxa, sqxr = sqrt(1 + xp^2), sqrt(1 + xa^2), sqrt(1 + xr^2) 
    # Analytical expression of (dr/du)(1/vr), that is always well-posed
    return (
        3xr * sqrt(
            (sqxr + sqxp) 
            * (sqxr + sqxa) 
            * (sqxp + sqxa)
            / (xr + xp) 
            / (xr + xa) 
            / (4 - u^2)
        )
        / (sqrt(2) * frequency_scale(model))
    )
end

function dΘdu(
    u::Float64,
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    bc = model.bc
    rp, ra = rpra_from_ae(a, e)
    r = radius_from_anomaly(u,a,e)
    #=

     the isochrone analytic Jacobian, Fouvry 21 G10
    =#
    xp, xa, xr = rp/bc, ra/bc, r/bc # Rescaled pericentre, apocentre, and radius
    sr = sqrt(1+xr^(2))
    
    Ω1, _ = frequencies_from_ae(a, e, model, params)
    sp,sa = spsa_from_ae(rp,ra,model)
    return (
        (3/sqrt(2))
        * (Ω1/frequency_scale(model))
        * (xr/sqrt(4-u^(2)))
        * (sqrt((sr+sp)*(sr+sa)*(sp+sa))/sqrt((xr+xp)*(xr+xa)))
    )
end


########################################################################
#
# (a, e) ↔ frequencies mappings
#
########################################################################
"""
for isochrone analytical version, see equations (G5-G7) in Fouvry&Prunet (2022)
"""
function αβ_from_ae(
    a::Float64,
    e::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    rp, ra = rpra_from_ae(a, e)
    xp, xa = (rp, ra) ./ model.bc
    sp, sa = spsa_from_rpra(rp, ra, model)
    return (2 / (sp + sa))^(3/2), (1 + (xp * xa) / ((1 + sp) * (1 + sa))) / 2
end

function ae_from_αβ(
    α::Float64,
    β::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    # (α,β) ↦ (E,L), analytic expressions specific to isochrone potential
    E, L = EL_from_αβ(α, β, model, params)
    # (E,L) ↦ (a,e), analytic expressions specific to isochrone potential
    a, e = ae_from_EL(E, L, model, params)
    return a, e
end

function ae_from_frequencies(
    Ω1::Float64,
    Ω2::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    α, β = αβ_from_frequencies(Ω1, Ω2, model)
    return ae_from_αβ(α, β, model, params)
end

function _β_from_α_circular(
    α::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)::Float64
    return 1/(1 + α^(2/3))
end

########################################################################
#
# (E,L) ↔ frequencies mappings
#
########################################################################
"""
for isochrone analytical version, see equations @TOCOMPLETE in Fouvry&Prunet (2022)
"""
function EL_from_αβ(
    α::Float64,
    β::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    return (
        0.5 * energy_scale(model) * α^(2/3), 
        momentum_scale(model) * (2β - 1) / (sqrt(β * (1 - β)))
    )
end


"""
analytic jacobian for the change of variables from E,L to α,β

@WARNING: given as a function of `α` and `β`, very confusing !!
"""
function EL_to_αβ_jacobian(
    α::Float64,
    β::Float64,
    model::AnalyticIsochrone,
    params::OrbitalParameters=OrbitalParameters()
)
    return abs(
        (1/6) 
        * energy_scale(model) 
        * momentum_scale(model) 
        / (
            α^(1/3) 
            * (β * (1 - β))^(3/2)
        )
    ) 
end



############################################################################################
#
#
# @IMPROVE: Function below not adapted to v2.0
#
#
############################################################################################
"""
# Returning the resonance frequency,
# along the circular orbit
"""
function _ωn(n1::Int64,n2::Int64,α::Float64)
    return n1*α + n2*α*_β_from_α_circular(α,model,params) # Output of the circular frequency
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