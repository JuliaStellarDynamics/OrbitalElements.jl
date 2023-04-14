"""
"""
function isochroneQROI(E::Float64,L::Float64,Ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochroneE0(bc,M,astronomicalG)
    Q = (E + L^(2)/(2.0*Ra^(2)))/scaleEnergy # Computing the Q parameter
    return Q # Output
end


"""
Saha distribution function
ra is the anisotropy radius
"""
function isochroneSahaDF(E::Float64,L::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    #scaleEnergy = isochroneE0(bc,M,astronomicalG)

    Q = isochroneQROI(E,L,ra,bc,M,astronomicalG)
    gamma = (bc/ra)^2

    prefactor = (M/((G*M*bc)^(3/2))) * (1/(128*sqrt(2)*pi)) * (sqrt(Q)/((1-Q)^4))

    term1 = 27 + 77gamma - (66+286gamma)*Q + (320+136gamma)*Q*Q -(240+32gamma)*Q*Q*Q + 64*Q*Q*Q*Q
    term2 = ((3*asin(sqrt(Q)))/sqrt(Q*(1-Q))) * ((-9 + 17gamma) + (28-44gamma)*Q + (16-8gamma)*Q*Q)

    return prefactor * (term1+term2)
end


"""
we also need dF/dQ and then we have the isochrone complete!
@ATTENTION, check the pi^3 factor here
@ATTENTION, check this answer against isochroneSahaDF
"""
function DFISOQROI(Q::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    gammaROI = (bc/ra)^2
    return (Mtot*sqrt(Q)*(27.0+77.0*gammaROI+2.0*Q*(-11.0*(3.0+13.0*gammaROI)+
    4.0*Q*(40+2.0*Q*(-15.0+4.0*Q-2.0*gammaROI)+17.0*gammaROI))+
    (3.0*(-9.0+4.0*Q*(7.0-2.0*Q*(-2.0+gammaROI)-11.0*gammaROI)+17.0*gammaROI)*
    asin(sqrt(Q)))/sqrt(-((-1.0+Q)*Q))))/(128.0*sqrt(2.0)*(bc*astronomicalG*M)^(3/2)*(pi)^(3)*(-1.0+Q)^(4))
end

"""
Saha distribution function derivative w.r.t. Q
ra is the anisotropy radius
"""
function isochroneSahadDFdQ(Q::Float64,ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    gammaROI = (bc/ra)^2
    return (M*((-1.0+Q)*(-128.0*gammaROI+Q*(-75.0+451.0*gammaROI+2.0*Q*
    (-659.0+387.0*gammaROI+4.0*Q*(90-23.0*gammaROI+2.0*Q*(-21.0+4.0*Q+2.0*gammaROI)))))-
    15.0*sqrt(-((-1.0+Q)*Q))*(5.0-13.0*gammaROI+4.0*Q*(-13.0+2.0*Q*(-2.0+gammaROI)+
    17.0*gammaROI))*asin(sqrt(Q))))/(256.0*sqrt(2.0)*(bc*astronomicalG*M)^(3/2)*(pi)^(3)*(-1.0+Q)^(6)*sqrt(Q))
end

"""
jacobian for converting dF/dQ dQ/dE -> dF/dE
"""
function isochronedQdEROI(E::Float64,L::Float64,Ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochroneE0(bc,M,astronomicalG)
    return 1.0/scaleEnergy # Output
end


"""
jacobian for converting dF/dQ dQ/dL -> dF/dL

@IMPROVE, does this need an extra factor of bc?
"""
function isochronedQdLROI(E::Float64,L::Float64,Ra::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    scaleEnergy = isochroneE0(bc,M,astronomicalG)
    return (L)/(scaleEnergy*Ra^(2)) # Output
end



"""isochroneisotropicDF(E[,bc,M,G])
the isotropic DF
"""
function isochroneisotropicDF(E::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)
    # get the relevant scale energy
    scaleEnergy = isochroneE0(bc,M,astronomicalG)

    # rescale dimensionless energy to goes between 0 and 1/2
    mE = E/scaleEnergy

    return (sqrt(mE)*M*(27.0+2.0*mE*(-1.0+4.0*mE)*(33.0+4.0*mE*(-7.0+2.0*mE))+
    (3.0*(-9.0+4.0*mE*(7.0+4.0*mE))*asin(sqrt(mE)))/sqrt(-((-1.0+mE)*mE))))/
    (128.0*sqrt(2.0)*(-1.0+mE)^(4)*(bc*astronomicalG*M)^(3/2)*(pi)^(3))

end


"""isochroneisotropicdDFdE(E[,bc,M,G])
the isotropic DF derivative w.r.t. E
"""
function isochroneisotropicdDFdE(E::Float64,bc::Float64=1.,M::Float64=1.,astronomicalG::Float64=1.)

    # get the relevant scale energy
    scaleEnergy = isochroneE0(bc,M,astronomicalG)

    # rescale dimensionless energy to goes between 0 and 1/2
    mE = E/scaleEnergy

    # return the isotropic derivative
    return (M*((-1.0+mE)*sqrt(mE)*(-75.0+2.0*mE*(-659.0+8.0*mE*(45.0+mE*(-21.0+4.0*mE))))+
    15.0*sqrt(1.0-mE)*(-5.0+4.0*mE*(13.0+4.0*mE))*asin(sqrt(mE))))/
    (256.0*sqrt(2.0)*scaleEnergy*(-1.0+mE)^(6)*(bc*astronomicalG*M)^(3/2)*(pi)^(3))

end
