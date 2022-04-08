"""
numerical inversion of (omega1,omega2) -> (a,e)
by brute-forcing the derivative steps domega1/da, domega1/de, deomega2/da, domega2/de

"""



"""ae_from_omega1omega2_brute

basic interpolation version to find (a,e) from (omega1,omega2) brute force derivatives.


"""
function ae_from_omega1omega2_brute(omega1::Float64,omega2::Float64,
                                    potential::Function,
                                    dpotential::Function,
                                    ddpotential::Function)
    #

    # get the circular orbit (maximum radius) for a given omega1,omega2. use the stronger constraint.
    aguess = Omega1circ_to_radius(omega1,dpotential,ddpotential)

    # then start from ecc=0.5 and take numerical derivatives
    eguess = 0.5
    f1,f2,df1da,df2da,df1de,df2de = compute_frequencies_ae_derivs(potential,dpotential,ddpotential,aguess,eguess)

    # 2d Newton Raphson inversion and find new steps

    return aguess,eguess
end
