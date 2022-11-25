"""Extremise.jl

helper routines for extremising functions (finding zeros)

"""

"""bisection(fun, xl, xu [, tolx, tolf])

A simple bisection algorithm, but it makes no allocations and is sufficiently fast: find the zero of a monotonic function.

"""
function bisection(fun::Function,
                   xl::Float64,
                   xu::Float64,
                   tolx::Float64=1000.0*eps(Float64),
                   tolf::Float64=1000.0*eps(Float64),
                   VERBOSE::Int64=0)
    """
    @ATTENTION, the tolerances are most likely overkill -- it may prevent convergence for high-order schemes

    """

    # order the input arguments
    if (xl > xu)
        xl, xu = xu, xl
    end

    # evaluate the function on the bracket
    fl, fu = fun(xl), fun(xu)

    if (abs(fl) <= tolf) # We have already found a solution on the left bracket
        return xl # Returning the left bracket
    end
    #####
    if (abs(fu) <= tolf) # We have already found a solution on the right bracket
        return xu # Returning the right bracket
    end
    #####
    @assert fl*fu < 0.0 "bisection: NOT A BRACKET"
    #####
    # Maximal number of steps to get to the tolerance in x
    nitermax = convert(Int64,ceil(log2(abs(xl-xu)/tolx)))
    for k = 1:nitermax # Bisection loop
        #####
        xm = (xl+xu)*0.5 # Middle value
        #####
        if (abs(xu-xl) <= tolx) # The considered bracket is smaller than the tolerance
            return xm # Returning the middle value
        end
        #####
        fm = fun(xm) # Value at the midpoint
        #####
        if (abs(fm) <= tolf) # The middle value is below the threshold
            return xm # Returning the middle value
        end
        #####
        # Otherwise, we iterate the bisection
        if (fm*fl < 0.0) # One root can be found between the left point and the middle point
            xu, fu = xm, fm # The upper point becomes the midpoint
        else
            xl, fl = xm, fm # The lower point becomes the midpoint
        end
        if k == nitermax
            if VERBOSE > 0
                println("OrbitalElements.Extremise.bisection: Maximal number of iteration reached.")
            end
            return xm
        end
    end

end


"""ExtremiseFunction(fun, xl, xu [,dx, tolf,VERBOSE])
Find the single extremum of a function (with monotonic derivative) between xl and xu
"""
function ExtremiseFunction(fun::Function,
                           xl::Float64=0.,
                           xu::Float64=1.,
                           dx::Float64=1.e-9,
                           tolf::Float64=1000.0*eps(Float64),
                           VERBOSE::Int64=0)

    # Finite differences derivative (only up +dx because xl can be 0. and fun not defined for x<0)
    dfun = x -> (fun(x+dx)-fun(x))/dx
    # Searching for derivative cancellation
    # Precision on x cannot be better than dx
    xm = try bisection(dfun,xl,xu,dx,tolf,VERBOSE) catch; (abs(fun(xu)) < abs(fun(xl))) ? xl : xu end
    return xm
end

"""ExtremiseFunctionNulCure(fun, xl, xu [, VERBOSE])
Find the extremum of a function between xl and xu.
Cure the possible nul (finite difference) derivative in xl or xu
"""
function ExtremiseFunctionNulCure(fun::Function,
                                  xl::Float64=0.,
                                  xu::Float64=1.,
                                  tolx::Float64=1000.0*eps(Float64),
                                  VERBOSE::Int64=0)

    dx = 1.e-3
    dfxl, dfxu = (fun(xl+dx)-fun(xl)), (fun(xl+dx)-fun(xl))

    while (dfxl == 0.) || (dfxu == 0.)
        dx *= 10.
        dfxl, dfxu = (fun(xl+dx)-fun(xl)), (fun(xl+dx)-fun(xl))
    end

    xm = ExtremiseFunction(fun,xl,xu,dx,0.,VERBOSE)

    if (xm == xl) || (xm == xu)
        return xm
    else
        xm = ExtremiseFunction(fun,xm-dx,xm+dx,tolx,0.,VERBOSE)
        return xm
    end
end
