"""
basic orbit transformations

"""

"""
    ae_from_rpra(rp, ra)

translate pericentre and apocentre to semi-major axis and eccentricity
"""
function ae_from_rpra(rp::Float64,ra::Float64)::Tuple{Float64,Float64}
    _checkdomain_rpra(rp, ra)
    return (rp + ra) / 2, (ra - rp) / (rp + ra)
end

"""
    rpra_from_ae(a, e)

translate semi-major axis and eccentricity to pericentre and apocentre
"""
function rpra_from_ae(a::Float64,e::Float64)::Tuple{Float64,Float64}
    _checkdomain_ae(a, e)
    return a * (1 - e), a * (1 + e)
end

"""
    rpra_to_ae_jacobian(rp,ra)

Jacobian of the (rp,ra) ↦ (a,e) mapping, i.e. |∂(a,e)/∂(rp,ra)|
"""
function rpra_to_ae_jacobian(rp::Float64,ra::Float64)::Float64
    _checkdomain_rpra(rp, ra)
    return 1 / (rp + ra)
end


"""
    ae_to_rpra_jacobian(a,e)

Jacobian of the (a,e) ↦ (rp,ra) mapping, i.e. |∂(rp,ra)/∂(a,e)|
"""
function ae_to_rpra_jacobian(a::Float64,e::Float64)::Float64
    _checkdomain_ae(a, e)
    return 2a
end