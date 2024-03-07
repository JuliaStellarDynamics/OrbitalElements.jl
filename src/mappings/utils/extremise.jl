"""Extremise.jl

helper routines for extremising functions (finding zeros)
"""

"""
    _bisection(fun, xl, xu[, tolx, tolf, VERBOSE])

A simple bisection algorithm, but it makes no allocations and is sufficiently fast:
find the zero of a monotonic function.

@ATTENTION: the tolerances are most likely overkill -- it may prevent convergence for
high-order schemes
"""
function _bisection(
    fun::Function,
    xl::Float64,
    xu::Float64;
    tolx::Float64=1000.0*eps(Float64),
    tolf::Float64=1000.0*eps(Float64)
)
    # order the input arguments
    if xl > xu
        xl, xu = xu, xl
    end
    # evaluate the function on the bracket
    fl, fu = fun(xl), fun(xu)
    #####
    if abs(fl) <= tolf # We have already found a solution on the left bracket
        return xl # Returning the left bracket
    end
    #####
    if abs(fu) <= tolf # We have already found a solution on the right bracket
        return xu # Returning the right bracket
    end
    #####
    @assert fl*fu < 0.0 "bisection: NOT A BRACKET"
    #####
    # Maximal number of steps to get to the tolerance in x
    nitermax = ceil(Int64, log2(abs(xl - xu) / tolx))
    for k = 1:nitermax # Bisection loop
        #####
        xm = (xl + xu) / 2 # Middle value
        #####
        if abs(xu-xl) <= tolx # The considered bracket is smaller than the tolerance
            return xm # Returning the middle value
        end
        #####
        fm = fun(xm) # Value at the midpoint
        #####
        if abs(fm) <= tolf # The middle value is below the threshold
            return xm # Returning the middle value
        end
        #####
        # Otherwise, we iterate the bisection
        if fm*fl < 0.0 # One root can be found between the left point and the middle point
            xu, fu = xm, fm # The upper point becomes the midpoint
        else
            xl, fl = xm, fm # The lower point becomes the midpoint
        end
        if k == nitermax
            warn("Maximal number of iterations reached in bisection.")
            return xm
        end
    end
end


"""
    _extremise(fun, xl, xu[, dx, tolf])

Find the extremum of a function (with monotonic derivative) between xl and xu.

@ASSUMPTION: `fun` is continuous and has a single extrema in `[xl,xu]`
"""
function _extremise(fun::Function,
    xl::Float64=0.,
    xu::Float64=1.;
    dx::Float64=1.e-9,
    tolf::Float64=1000.0*eps(Float64)
)
    # Finite differences derivative (only up +dx because xl can be 0. and 
    # fun not defined for x<0)
    dfun = x -> (fun(x + dx) - fun(x)) / dx
    # Searching for derivative cancellation
    # Precision on x cannot be better than dx
    xm = try 
        _bisection(dfun, xl, xu; tolx=dx, tolf=tolf) 
    catch; 
        abs(fun(xu)) < abs(fun(xl)) ? xl : xu 
    end
    return xm
end

"""
    _extremise_noedges(fun,xl,xu[,tolx,VERBOSE])

as [`_extremise`](@ref) but cure the possible nul derivative in xl or xu.
"""
function _extremise_noedges(
    fun::Function,
    xl::Float64=0.,
    xu::Float64=1.;
    tolx::Float64=1000.0*eps(Float64)
)
    dx = 1e-3
    dfxl, dfxu = fun(xl + dx) - fun(xl), fun(xl + dx) - fun(xl)

    while (dfxl == 0.) || (dfxu == 0.)
        dx *= 10
        dfxl, dfxu = fun(xl + dx) - fun(xl), fun(xl + dx) - fun(xl)
    end

    xm = _extremise(fun, xl, xu, dx=dx, tolf=0.)

    if xm == xl || xm == xu
        return xm
    else
        xm = _extremise(fun, xm-dx, xm+dx, dx=tolx, tolf=0.)
        return xm
    end
end
