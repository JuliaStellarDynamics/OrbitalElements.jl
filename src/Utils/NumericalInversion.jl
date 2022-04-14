"""
numerical inversion of (omega1,omega2) -> (a,e)
by brute-forcing the derivative steps domega1/da, domega1/de, deomega2/da, domega2/de

"""



"""ae_from_omega1omega2_brute

basic interpolation version to find (a,e) from (omega1,omega2) brute force derivatives.

@IMPROVE add escape for circular orbits
"""
function ae_from_omega1omega2_brute(omega1::Float64,omega2::Float64,
                                    potential::Function,
                                    dpotential::Function,
                                    ddpotential::Function,
                                    eps::Float64=1*10^(-6),
                                    maxiter::Int64=10000)
    #

    # get the circular orbit (maximum radius) for a given omega1,omega2. use the stronger constraint.
    acirc = Omega1circ_to_radius(omega1,dpotential,ddpotential)


    # check to make sure we aren't very close to circular (radial isn't a problem)
    f1circ,f2circ = compute_frequencies_ae(potential,dpotential,ddpotential,acirc,0.0)
    #if (((omega1 - f1circ)^2 + (omega2 - f2circ)^2) < eps^2)
    #    return acirc,0.0
    #end


    # then start from ecc=0.5 and take numerical derivatives
    aguess = acirc
    eguess = 0.5
    f1,f2 = compute_frequencies_ae(potential,dpotential,ddpotential,aguess,eguess)

    # 2d Newton Raphson inversion and find new steps
    iter = 0
    while (((omega1 - f1)^2 + (omega2 - f2)^2) > eps^2)

        f1,f2,df1da,df2da,df1de,df2de = OrbitalElements.compute_frequencies_ae_derivs(potential,dpotential,ddpotential,aguess,eguess)

        jacobian = [df1da df1de ; df2da df2de]
        step = jacobian \ (-([f1;f2] - [omega1 ; omega2]))

        aguess,eguess = aguess + step[1],eguess + step[2]

        iter += 1
        if iter > maxiter
            break
        end
    end

    return aguess,eguess
end
