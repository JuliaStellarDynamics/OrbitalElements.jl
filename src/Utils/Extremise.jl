"""Extremise.jl


"""

"""ExtremiseFunction

Find the single extremum of a function between minval and maxval

Accuracy will be set by deps, 1/2^neps:
  can never do better than evaluating the function on a grid this fine.

Requires exactly 2*neps evaluations of the function
"""
function ExtremiseFunction_legacy(func::Function,neps::Int64=32,minval::Float64=0.,maxval::Float64=1.;verbose::Bool=false,fullreturn::Bool=false)

    if neps > 50
        # this should be able to do 53 for double precision, but I've found some bugs before that.
        neps = 50
        println("OrbitalElements.Extremise.jl: ExtremiseFunction- reached maximum precision.")
    end

    # set the finite difference derivative step size as the requested minimum
    deps = 1/2^(neps)

    width = maxval - minval

    # check the endpoint derivatives
    left_end_derivative  = func(minval+deps)-func(minval)
    right_end_derivative = func(maxval) - func(maxval-deps)

    if left_end_derivative*right_end_derivative > 0
        if verbose
            println("OrbitalElements.Extremise.jl: Input function is monotonic")
        end
        return -1
    end

    left_derivative  = left_end_derivative
    right_derivative = right_end_derivative
    leftmin  = minval
    rightmax = minval+width

    # set flag for final precision step counter
    rflag = false

    for iter=1:neps
        # assign a new midpoint
        midpoint = leftmin + width/(2^iter)
        midfunc = func(midpoint)

        # compute finite difference derivatives on either side of the midpoint
        mid_left_derivative  = midfunc - func(midpoint-deps)
        mid_right_derivative = func(midpoint+deps) - midfunc

        # check here for NaN values:

        # consider a special case of equal derivatives to block extra evaluations
        #print(mid_left_derivative,'=',mid_right_derivative)
        # if first step, turn up deps and check that we aren't in a weird locality
        if mid_left_derivative == -mid_right_derivative
            # we are directly on the extremum!
            if iter==1
                if verbose
                    println("OrbitalElements.Extremise.jl: Turning up deps")
                    deps = 1/2^(floor(Int, neps/2))
                    mid_left_derivative  = midfunc - func(midpoint-deps)
                    mid_right_derivative = func(midpoint+deps) - midfunc

                    # what if we are directly at the midpoint?
                    if mid_left_derivative == -mid_right_derivative
                        if verbose
                            println("OrbitalElements.Extremise.jl: Directly centred at $midpoint")
                        end
                        return midpoint
                    end

                end
                continue
            else
                if verbose
                    println("OrbitalElements.Extremise.jl: Directly centred at $midpoint")
                end
                return midpoint
            end
        end

        # check if on left side of the interval
        if mid_left_derivative*left_derivative < 0
            # set the right side (maximum) to be the midpoint
            rightmax = midpoint
            # set the right derivative to be the old midpoint left derivative
            right_derivative = mid_left_derivative
        end

        # check if on the right side of the interval
        # if not on left, must be on the right side!
        if mid_right_derivative*right_derivative < 0
            # set the left side (minimum) to be the midpoint
            leftmin = midpoint
            # set the left derivative to be the old midpoint right derivative
            left_derivative = mid_right_derivative

            # check if we are at the final step
            if (iter == neps)
                rflag = true
            end
        end

        # watch the convergence fly by...
        if verbose
            println("OrbitalElements.Extremise.jl: #$iter: best=$midpoint, min=$leftmin, max=$rightmax")
        end

        # can add some fancy check here for oscillatory convergence, if we want

        # otherwise, return the value at the maximum number of iterations
        # (with an advancement to the next midpoint)
        if iter == neps

            #check to make sure we aren't at the boundaries
            #print(midpoint,' ',maxval-deps,'\n')
            if (midpoint>=(maxval-deps)) | (midpoint<=(minval+deps))
                return -1
            end

            # take one last midpoint based on which side the derivative is on
            if rflag
                if fullreturn
                    return midpoint + 1/(2^(neps+1)),leftmin,rightmax
                else
                    return midpoint + 1/(2^(neps+1))
                end
            else
                if fullreturn
                    return midpoint - 1/(2^(neps+1)),leftmin,rightmax
                else
                    return midpoint - 1/(2^(neps+1))
                end
            end

        end
    end
end


"""bisection(fun, xl, xu [, tolx, tolf])
used to find zero of a monotonic function.
A simple bisection algorithm, but it makes no allocations and is sufficiently fast
@ATTENTION, the tolerances are most likely overkill -- it may prevent convergence for high-order schemes
"""
function bisection(fun::Function,
                   xl::Float64,
                   xu::Float64;
                   tolx::Float64=1000.0*eps(Float64),
                   tolf::Float64=1000.0*eps(Float64),
                   nitermax::Int64=100)

    if (xl > xu)
        xl, xu = xu, xl # Ordering the input arguments
    end
    #####
    fl, fu = fun(xl), fun(xu) # Evaluating the function on the bracket
    #####
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
            println("Warning : Maximal number of iteration reached in bisection.")
            return xm
        end
    end

end


"""ExtremiseFunction
accepts a function to find the zero, or to maximise the derivative
"""
function ExtremiseFunction(fun::Function,
                           xl::Float64=0.,
                           xu::Float64=1.,
                           dx::Float64=1e-8;
                           tolx::Float64=1000.0*eps(Float64),
                           tolf::Float64=1000.0*eps(Float64))

    dfun = x -> (fun(x+dx)-fun(x))/(dx)
    #xext = try bisection(dfun,xl,xu;tolx=tolx,tolf=tolf); catch -1
    xext = bisection(dfun,xl,xu;tolx=tolx,tolf=tolf)

    return xext
end
