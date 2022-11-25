

"""
function to perform a Simpson's 1/3 composite integration between -1<u<1
@ATTENTION, the number of steps K must be even
fun is a function that takes u as an argument and returns the value of some function at that point
this function can also take a function that passes an array and push them all forward at once.
"""
function UnitarySimpsonIntegration(fun::Function,K::Int64)

    # specify the step size
    δu = 2.0 / K

    # set up the launching point
    u   = -1.0
    res = δu .* fun(u) ./ 6.0

    # add contributions from the interior points
    for k = 1:(K-1)
        u   += 0.5 * δu
        res = res .+ 2.0 .* δu .* fun(u) ./ 3.0
        u   += 0.5 * δu
        res = res .+ 1.0 .* δu .* fun(u) ./ 3.0
    end

    # add the final points
    u   += 0.5 * δu
    res = res .+ 2.0 .* δu .* fun(u) ./ 3.0
    u   += 0.5 * δu
    res = res .+ 1.0 .* δu .* fun(u) ./ 6.0

    # return the value(s)
    return res
end
