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
    NumericalDerivativePoints(a,e,da,de,tola,tole)
    
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
"""
function NumericalDerivativePoints(a::Float64,e::Float64,
                                    da::Float64,de::Float64,
                                    tola::Float64,tole::Float64)::Tuple{Float64,Float64,Float64,Float64}

    # Usual points
    ap = a + da
    ep = e + de
    # Check for borders and tolerance limits and adapt 
    if (a < tola) && (ap >= tola)
        if a - da <= 0.
            println("WARNING: Too low tolerance on semimajor axis compared to the numerical derivative step.")
        else
            ap = a - da
            da *= -1.0
        end
    end
    if (e < tole) && (ep >= tole)
        if e - de < 0.
            println("WARNING: Too low tolerance on eccentricity compared to the numerical derivative step.")
        else
            ep = e - de
            de *= -1.0
        end
    elseif ep > 1.0 - tole
        if (ep > 1.0) || (e < 1.0 - tole)
            if (ep > 1.0) && (e - de < 1.0 - tole)
                println("WARNING: Too low tolerance on eccentricity compared to the numerical derivative step.")
            end
            ep = e - de
            de *= -1.0
        end
    end

    return ap, da, ep, de
end


"""
    NumericalDerivativeAE(fun,a,e,params)
    
compute the numerical derivative of any function of (a,e).
"""
function NumericalDerivativeAE(fun::F0,
                               a::Float64,e::Float64,
                               params::OrbitalParameters=OrbitalParameters()) where {F0 <: Function}
    
    # Numerical derivative points
    tola, tole = params.TOLA, EccentricityTolerance(a,params.rc,params.TOLECC)
    ap, da, ep, de = NumericalDerivativePoints(a,e,params.da,params.de,tola,tole)

    # Derivation outside the integral
    floc = fun(a,e)

    # For a derivatives
    fap = fun(ap,e)
    
    # For e derivatives
    fep = fun(a,ep)

    ∂f∂a = (fap.-floc)./da
    ∂f∂e = (fep.-floc)./de

    return floc, ∂f∂a, ∂f∂e
end

########################################################################
#
# Interpolation
#
########################################################################

function Interpolation1stOrder(x::Float64,
                               x0::Float64,y0::Float64,
                               x1::Float64,y1::Float64)::Float64

    return (x*y0 - x1*y0 - x*y1 + x0*y1)/(x0 - x1)
end

function Interpolation2ndOrder(x::Float64,
                               x0::Float64,y0,
                               x1::Float64,y1,
                               x2::Float64,y2)

    return ((x - x2).*((x - x1)*(x1 - x2).*y0 .+ (x - x0)*(-x0 + x2).*y1) .+ (x - x0)*(x - x1)*(x0 - x1).*y2)./((x0 - x1)*(x0 - x2)*(x1 - x2))
end