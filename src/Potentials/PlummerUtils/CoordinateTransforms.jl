


#
#
# transformations between coordinate definitions
#

"""
translate From (sp,sa) to (rp,ra)
"""
function SpSaFromRpRa(rp::Float64,ra::Float64,bc::Float64)
    sp,sa = sqrt(1+(rp/bc)^2),sqrt(1+(ra/bc)^2)
    return sp,sa
end

"""
translate From (rp,ra) to (sp,sa)
"""
function RpRaFromSpSa(sp::Float64, sa::Float64, bc::Float64)
    rp = bc * sqrt(abs(sp^2 - 1.0))
    ra = bc * sqrt(abs(sa^2 - 1.0))
    return  rp, ra
end


"""
translate From (u,a,e) to s
"""
function SFromUAE(u::Float64, sma::Float64, ecc::Float64)
    fu = u * (1.5 - 0.5*u^2)
    return sma * (1.0 + ecc * fu)
end

"""
translate From s to r
"""
function RFromS(s::Float64,bc::Float64)
    return bc * sqrt(abs(s^2 - 1.0))
end

"""
translate From (u,a,e) to r
"""
function RFromUAE(u::Float64, a::Float64, e::Float64, bc::Float64)
    return RFromS(SFromUAE(u,a,e),bc)
end

"""
translate From (u,rp,ra) to r
"""
function RFromURpRa(u::Float64, rp::Float64, ra::Float64, bc::Float64)
    sp,sa = SpSaFromRpRa(rp,ra,bc)

    # compute a,e for the anomaly: equation F8 of Tep+ (2022)
    a,e = AEFromRpRa(sp,sa)
    return RFromS(SFromUAE(u,a,e),bc)
end


function AEFromSpSa(sp::Float64, sa::Float64)
    sma = (sa + sp)/2
    ecc = (sa - sp)/(sa + sp)
    return sma, ecc
end
