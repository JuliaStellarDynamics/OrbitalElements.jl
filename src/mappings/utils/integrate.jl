

"""
    _integrate_simpson(fun, K)

perform a Simpson's 1/3 composite integration over [-1, 1]

`fun` is a function of one argument but can return multiples values.
"""
function _integrate_simpson(fun::Function, K::Int64)

    # specify the step size
    δu = 2.0 / K
    # set up the launching point
    u   = -1.0
    res = δu .* fun(u) ./ 6.0

    # add contributions from the interior points
    for k in 1:(K-1)
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

    return res
end
