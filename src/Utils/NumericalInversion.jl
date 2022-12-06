"""
numerical inversion of (Ω₁,Ω₂) -> (a,e)
by brute-forcing the derivative increments dΩ₁/da, dΩ₁/de, deΩ₂/da, dΩ₂/de

VERBOSE rules:
0: no printing
2: failure cases diagnostics
3: iteration counters

"""

function inverse2Dlinear(a::Float64,b::Float64,
                         c::Float64,d::Float64,
                         y1::Float64,y2::Float64)::Tuple{Float64,Float64}

    deta = a*d - b*c
    if deta == 0.
        return 0., 0.
    else
        return (d*y1 - b*y2)/deta, (a*y2 - c*y1)/deta
    end
end



"""AEFromΩ1Ω2Brute(Ω₁,Ω₂,ψ,dψ,d2ψ,d3ψ[,eps,maxiter,TOLECC,TOLA,da,de,VERBOSE])

basic Newton-Raphson algorithm to find (a,e) from (Ω₁,Ω₂) brute force derivatives.

"""
@inline function AEFromΩ1Ω2Brute(Ω₁::Float64,Ω₂::Float64,
                         ψ::Function,
                         dψ::Function,
                         d2ψ::Function,
                         d3ψ::Function,
                         d4ψ::Function,
                         params::OrbitsParameters)::Tuple{Float64,Float64,Int64,Float64}
    """
    @IMPROVE add escape for circular orbits

    """
    # get the circular orbit (maximum radius) for a given Ω₁,Ω₂. use the stronger constraint.
    acirc = RcircFromΩ1circ(Ω₁,dψ,d2ψ,params.rmin,params.rmax)

    # then start from ecc=0.5 and take numerical derivatives
    aguess = acirc
    eguess = 0.5

    f1,f2,df1da,df2da,df1de,df2de = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,aguess,eguess,params)

    tol = (Ω₁ - f1)^2 + (Ω₂ - f2)^2
    if (tol < (params.invε)^2)
        return aguess, eguess, iter, tol
    end

    # 2d Newton Raphson inversion and find new increments
    for iter = 1:params.ITERMAX

        increment1, increment2 = inverse2Dlinear(df1da,df1de,df2da,df2de,Ω₁-f1,Ω₂-f2) 

        # If non inversible 
        if (increment1 == 0.) && (increment2 == 0.)
            break
        end
        
        # If end point not in the domain ( a >= 0, e in [0,1] )
        if ((eguess == 0.) && (increment2 < 0.)) || ((eguess == 1.) && (increment2 > 0.))
            increment2 = 0.
        end
        while (aguess + increment1 < 0.) || (eguess + increment2 < 0.) || (eguess + increment2 > 1.) 
            increment1 /= 2
            increment2 /= 2
        end

        # Update guesses
        aguess += increment1
        eguess += increment2

        f1,f2,df1da,df2da,df1de,df2de = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,aguess,eguess,params)

        tol = (Ω₁ - f1)^2 + (Ω₂ - f2)^2
        if (tol < (params.invε)^2)
            return aguess, eguess, iter, tol
        end
    end

    return aguess,eguess,params.ITERMAX,tol
end


"""AEFromJLBrute(Ω₁,Ω₂,ψ,dψ,d2ψ,d3ψ[,eps,maxiter,TOLECC,TOLA,da,de,VERBOSE])

basic Newton-Raphson algorithm to find (a,e) from (Jᵣ,L) brute force derivatives.

"""
@inline function AEFromJLBrute(J::Float64,L::Float64,
                         ψ::Function,
                         dψ::Function,
                         d2ψ::Function,
                         d3ψ::Function,
                         d4ψ::Function,
                         params::OrbitsParameters)::Tuple{Float64,Float64,Int64,Float64}

    # get the circular orbit (maximum radius) for a given angular momentum.
    acirc = RcircFromL(L,dψ,params.rmin,params.rmax)

    # then start from ecc=0.5 and take numerical derivatives
    aguess = acirc
    eguess = 0.5

    Jguess, Lguess, dJgda, dLgda, dJgde, dLgde = ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,d3ψ,aguess,eguess,params)

    tol = (Jguess - J)^2 + (Lguess - L)^2
    if (tol < (params.invε)^2)
        return aguess, eguess, iter, tol
    end

    # 2d Newton Raphson inversion and find new increments
    for iter = 1:params.ITERMAX

        increment1, increment2 = inverse2Dlinear(dJgda,dJgde,dLgda,dLgde,J-Jguess,L-Lguess)

        # If non inversible 
        if (increment1 == 0.) && (increment2 == 0.)
            break
        end
        
        # If end point not in the domain ( a >= 0, e in [0,1] )
        if ((eguess == 0.) && (increment2 < 0.)) || ((eguess == 1.) && (increment2 > 0.))
            increment2 = 0.
        end
        while (aguess + increment1 < 0.) || (eguess + increment2 < 0.) || (eguess + increment2 > 1.) 
            increment1 /= 2
            increment2 /= 2
        end

        # Update guesses
        aguess += increment1
        eguess += increment2

        Jguess, Lguess, dJgda, dLgda, dJgde, dLgde = ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,d3ψ,aguess,eguess,params)

        tol = (Jguess - J)^2 + (Lguess - L)^2
        if (tol < (params.invε)^2)
            return aguess, eguess, iter, tol
        end
    end

    return aguess,eguess,params.ITERMAX,tol
end