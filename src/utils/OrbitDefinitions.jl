"""
basic orbit transformations

"""

"""
    AEfromRpRa(rp,ra)

translate pericentre and apocentre to semi-major axis and eccentricity
"""
function AEFromRpRa(rp::Float64,ra::Float64)::Tuple{Float64,Float64}

    return (rp+ra)/2, (ra-rp)/(rp+ra)
end

"""
    RpRafromAE(a,e)

translate semi-major axis and eccentricity to pericentre and apocentre
"""
function RpRaFromAE(a::Float64,e::Float64)::Tuple{Float64,Float64}

    return a*(1-e), a*(1+e)
end

"""
    JacRpRaToAE(rp,ra)

Jacobian of the (rp,ra) ↦ (a,e) mapping, i.e. |∂(a,e)/∂(rp,ra)|
"""
function JacRpRaToAE(rp::Float64,ra::Float64)::Float64

    return 1.0 / (rp + ra)
end


"""
    JacAEToRpRa(a,e)

Jacobian of the (a,e) ↦ (rp,ra) mapping, i.e. |∂(rp,ra)/∂(a,e)|
"""
function JacAEToRpRa(a::Float64,e::Float64)::Float64

    return 2a
end