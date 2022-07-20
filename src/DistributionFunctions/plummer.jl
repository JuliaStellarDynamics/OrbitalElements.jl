

"""
Q(E,L), the variable for anisotropy distribution functions.

when ra->infty, this is isotropic.
"""
function plummer_Q_ROI(E::Float64,L::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = - astronomicalG * M / bc
    Q = (E + L^(2)/(2.0*ra^(2)))/scaleEnergy # Computing the Q parameter
    return Q # Output
end

"""
jacobian for converting dF/dQ dQ/dE -> dF/dE
"""
function plummer_dQdE_ROI(E::Float64,L::Float64,Ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = - astronomicalG * M / bc
    return 1.0/scaleEnergy # Output
end

"""
jacobian for converting dF/dQ dQ/dL -> dF/dL
"""
function plummer_dQdL_ROI(E::Float64,L::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = - astronomicalG * M / bc
    return (L)/(scaleEnergy*ra^(2)) # Output
end

function plummer_dFdQ(Q::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    return sqrt(-Q)
end
