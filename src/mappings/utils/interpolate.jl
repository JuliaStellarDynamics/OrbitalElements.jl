########################################################################
#
# Interpolation
#
########################################################################

function _interpolate_order_1(
    x::Float64,
    x0::Float64,
    y0,
    x1::Float64,
    y1
)
    return ((x - x1) .* y0 - (x - x0) .* y1)/(x0 - x1)
end

function _interpolate_order_2(
    x::Float64,
    x0::Float64,
    y0,
    x1::Float64,
    y1,
    x2::Float64,
    y2
)
    return (
        (
            (x - x2) .* ((x - x1) * (x1 - x2) .* y0 .+ (x - x0) * (-x0 + x2) .* y1) 
            .+ (x - x0) * (x - x1) * (x0 - x1) .* y2
        )
        ./ (
            (x0 - x1) * (x0 - x2) * (x1 - x2)
        )
    )
end