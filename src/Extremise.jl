"""Extremise.jl

This should almost certainly be it's own julia tool, but we can keep it here for now.

"""

"""extremise_function

Find the single extremum of a function between minval and maxval

Accuracy will be set by deps, 1/2^neps:
  can never do better than evaluating the function on a grid this fine.

Requires exactly 2*neps evaluations of the function


"""
function extremise_function(func::Function,neps::Int64=32,minval::Float64=0.,maxval::Float64=1.,verbose::Bool=false)

    if neps > 50
        # this should be able to do 53 for double precision, but I've found some bugs before that.
        neps = 50
        print("extremise_function- reached maximum precision.\n")
    end

    deps = 1/2^(neps)

    width = maxval - minval

    # check the endpoint derivatives
    left_end_derivative = func(minval+deps)-func(minval)
    right_end_derivative = func(maxval) - func(maxval-deps)

    if left_end_derivative*right_end_derivative > 0
        if verbose
            print("Monotonic\n")
        end
        return -1
    end

    left_derivative  = left_end_derivative
    right_derivative = right_end_derivative
    leftmin = 0.
    rightmax = width

    # set flag for final precision step counter
    rflag = false

    for iter=1:neps
        # assign a new midpoint
        midpoint = leftmin + width/(2^iter)
        midfunc = func(midpoint)

        # compute finite difference derivatives on either side of the midpoint
        mid_left_derivative  = midfunc - func(midpoint-deps)
        mid_right_derivative = func(midpoint+deps) - midfunc

        # consider a special case of equal derivatives to block extra evaluations
        #print(mid_left_derivative,'=',mid_right_derivative)
        # if first step, turn up deps and check that we aren't in a weird locality
        if mid_left_derivative == -mid_right_derivative
            # we are directly on the extremum!
            if iter==1
                if verbose
                    print("Turning up deps\n")
                    deps = 1/2^(floor(Int, neps/2))
                    mid_left_derivative  = midfunc - func(midpoint-deps)
                    mid_right_derivative = func(midpoint+deps) - midfunc

                    # what if we are directly at the midpoint?
                    if mid_left_derivative == -mid_right_derivative
                        if verbose
                            print("Directly centred at ",midpoint,"\n")
                        end
                        return midpoint
                    end

                end
                continue
            else
                if verbose
                    print("Directly centred at ",midpoint,"\n")
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
        #else
        if mid_right_derivative*right_derivative < 0 # if not on left, must be on the right side
            # set the left side (minimum) to be the midpoint
            leftmin = midpoint

            left_derivative = mid_right_derivative
            if (iter == neps)
                rflag = true
            end
        end

        # watch the convergence fly by...
        if verbose
            print(iter,' ',midpoint,' ',leftmin,' ',rightmax,'\n')
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
                return midpoint + 1/(2^(neps+1))
            else
                return midpoint - 1/(2^(neps+1))
            end

        end
    end
end
