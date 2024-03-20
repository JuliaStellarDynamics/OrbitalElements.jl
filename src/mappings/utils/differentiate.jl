"""
DerivationInterpolation.jl

Useful functions for numerical derivation or interpolations

"""


########################################################################
#
# Derivation
#
########################################################################

"""
    _derivatives_points_ae(a,e,da,de,tola,tole)

Points to use for numerical derivative w.r.t a and e
depending on the location (switch close to border and close to cut off)
default mode being 1st order right derivative [x,x+dx].

Points are structured as follow:
 [+da]
(ap, e)
   ↑
(a , e) → (a ,ep)
           [+de]

Output order :
semimajor axis derivative info followed by eccentricity
  da      de
ap, da, ep, de
(da and de could be switched to negative values)

@WARNING: Important assumption here
    → tolerances switch are made with exclusive lower or greater boundary conditions
    i.e. switch points are part of the standard case (not the border ones)

@IMPROVE: Can we do something about the warnings?
@IMPROVE: Should be a function of the dimensionless semimajor axes
"""
function _derivatives_points_ae(
    a::Float64,
    e::Float64,
    da::Float64,
    de::Float64,
    tola::Float64,
    tole::Float64
)::Tuple{Float64,Float64,Float64,Float64}
    # Usual points
    ap = a + da
    ep = e + de
    # Check for borders and tolerance limits and adapt
    if a < tola && ap >= tola
        if a - da <= 0.
            println(
                "Too low tolerance on semimajor axis for the numerical derivative step."
            )
        else
            ap = a - da
            da *= -1.0
        end
    end
    if e < tole && ep >= tole
        if e - de < 0.
            println(
                "Too low tolerance on eccentricity for the numerical derivative step."
            )
        else
            ep = e - de
            de *= -1.0
        end
    elseif ep > 1.0 - tole
        if ep > 1.0 || (e < 1.0 - tole)
            if ep > 1.0 && e - de < 1.0 - tole
                println(
                    "Too low tolerance on eccentricity for the numerical derivative step."
                )
            else
                ep = e - de
                de *= -1.0
            end
        end
    end

    return ap, da, ep, de
end


"""
    _derivatives_ae(fun, a, e, params)

compute the numerical derivative of any function of semimajor axis and eccentricity using 
naive right finite differences.

@IMPROVE: diffentiation scheme is hard-coded. Make it a parameter
@IMPROVE: Should be a function of the dimensionless semimajor axis
"""
function _derivatives_ae(
    fun::F0,
    a::Float64,
    e::Float64,
    params::OrbitalParameters=OrbitalParameters()
) where {F0 <: Function}
    # Numerical derivative points
    # @IMPROVE: For now patched with the parameter rc, but the semimajor axes should be 
    # adimensional
    tola, tole = params.TOLA,_eccentricity_tolerance(a/params.rc, params.TOLECC)
    ap, da, ep, de = _derivatives_points_ae(a, e, params.da, params.de, tola, tole)
    
    floc = fun(a, e) # Both derivatives
    fap = fun(ap, e) # For a derivatives
    fep = fun(a, ep) # For e derivatives

    ∂f∂a = (fap .- floc) ./ da
    ∂f∂e = (fep .- floc) ./ de

    return floc, ∂f∂a, ∂f∂e
end