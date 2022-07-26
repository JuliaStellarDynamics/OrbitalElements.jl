

"""
Q(E,L), the variable for anisotropy distribution functions.

when ra->infty, this is isotropic.
"""
function plummer_ROI_Q(E::Float64,L::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = - astronomicalG * M / bc
    Q = (E + (L^(2)*(bc^2))/(2.0*ra^(2)))/scaleEnergy # Computing the Q parameter
    return Q # Output
end

"""
jacobian for converting dF/dQ dQ/dE -> dF/dE
"""
function plummer_ROI_dQdE(E::Float64,L::Float64,Ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = - astronomicalG * M / bc
    return 1.0/scaleEnergy # Output
end

"""
jacobian for converting dF/dQ dQ/dL -> dF/dL
"""
function plummer_ROI_dQdL(E::Float64,L::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = - astronomicalG * M / bc
    return (L)/(scaleEnergy*ra^(2)) # Output
end

"""
the anisotropic distribution function from PlummerPlus
"""
function plummer_ROI_DF(E::Float64,L::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    Q = plummer_ROI_Q(E,L,ra,bc,M,astronomicalG)

    # dimensionless anisotropy radius
    gamma = (bc/ra)^2

    prefactor = (M/((astronomicalG*M*bc)^(3/2))) * (48*sqrt(3)/(7*(pi^3)))

    return prefactor * (Q^(7/2)) * (1-gamma+7gamma/(16*(Q^2)))

end


"""
the derivative of the anisotropic distribution function
"""
function plummer_ROI_dFdQ(Q::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    #Q = plummer_ROI_Q(E,L,ra,bc,M,astronomicalG)

    # dimensionless anisotropy radius
    gamma = (bc/ra)^2

    prefactor = (M/((astronomicalG*M*bc)^(3/2))) * (48*sqrt(3)/(7*(pi^3)))

    return prefactor * sqrt(Q) * ( ((7/2)-(7gamma/2))*Q*Q + (21gamma/32))

end


"""
the isotropic distribution function
"""
function plummer_ISO_DF(E::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    scaleEnergy = - astronomicalG * M / bc

    prefactor = (M/((astronomicalG*M*bc)^(3/2))) * (48*sqrt(3)/(7*(pi^3)))

    return prefactor * (E/scaleEnergy)^(7/2)

end


"""
the derivative of the isotropic distribution function
"""
function plummer_ISO_dFdE(E::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    scaleEnergy = - astronomicalG * M / bc

    prefactor = (M/((astronomicalG*M*bc)^(3/2))) * (24*sqrt(3)/(7*(pi^3)))

    return prefactor * (E/scaleEnergy)^(5/2)

end



"""
calculate n.dF/dE, for the isotropic plummer case.
"""
function plummer_ISO_ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    dFdE = plummer_ISO_dFdE(E,bc,M,astronomicalG)

    # Value of n.dF/dJ
    res = ndotOmega*dFdE

    return res
end

"""
calculate n.dF/dE, for the osipkov-merritt plummer case.

@IMPROVE, make a distribution function fix. perhaps only integrate in regions where the DF is positive?
"""
function plummer_ROI_ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.,ra::Float64=1.)

    Q = plummer_ROI_Q(E,L,Ra,bc,M,astronomicalG)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    # Value of dF/dQ
    dFdQ = plummer_ROI_dFdQ(Q,ra,bc,M,astronomicalG)

    # Values of dQ/dE, dQ/dL
    dQdE, dQdL = plummer_ROI_dQdE(E,L,ra,bc,M,astronomicalG), plummer_ROI_dQdL(E,L,ra,bc,M,astronomicalG)

    # Value of n.dF/dJ
    res = dFdQ*(dQdE*ndotOmega + n2*dQdL)

    return res

end
