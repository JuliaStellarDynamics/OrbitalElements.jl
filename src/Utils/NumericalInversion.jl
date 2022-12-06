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
    da, de = params.da, params.de
    TOLECC = params.TOLECC
    # get the circular orbit (maximum radius) for a given Ω₁,Ω₂. use the stronger constraint.
    acirc = RcircFromΩ1circ(Ω₁,dψ,d2ψ,params.rmin,params.rmax)

    # then start from ecc=0.5 and take numerical derivatives
    aguess = acirc
    eguess = 0.5
    f1,f2 = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,aguess,eguess,params)

    # 2d Newton Raphson inversion and find new increments
    iter = 0
    while (((Ω₁ - f1)^2 + (Ω₂ - f2)^2) > (params.invε)^2)

        f1,f2,df1da,df2da,df1de,df2de = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,aguess,eguess,params)

        # one break: negative frequencies when getting very close to the centre.
        # define the failure mode: return circular orbit at the minimum size
        if (f1 < 0.0) | (f2 < 0.0)
            return da,0.0,iter+1,1.
        end

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

        iter += 1

        if iter > params.ITERMAX
            break
        end

    end


    #(VERBOSE > 2) && println("OrbitalElements.NumericalInversion.AEFromΩ1Ω2Brute: niter=",iter)

    finaltol = ((Ω₁ - f1)^2 + (Ω₂ - f2)^2)

    # check here to not allow bad values?
    if isnan(aguess) | isnan(eguess)
        #(VERBOSE > 1) && println("OrbitalElements.NumericalInversion.AEFromΩ1Ω2Brute: failed for inputs (Ω₁,Ω₂)=($Ω₁,$Ω₂).")

        # return a semi-equivalent circular orbit, as the failure mode is mostly very small orbits
        return acirc,0.0,-1,finaltol

    else
        return aguess,eguess,iter,finaltol
    end
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
    """
    @IMPROVE add escape for circular orbits

    """
    da, de = params.da, params.de
    TOLECC = params.TOLECC
    # get the circular orbit (maximum radius) for a given Ω₁,Ω₂. use the stronger constraint.
    acirc = RcircFromL(L,dψ,params.rmin,params.rmax)

    # then start from ecc=0.5 and take numerical derivatives
    aguess = acirc
    eguess = 0.5

    Jguess,Lguess = ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,aguess,eguess,params)

    # 2d Newton Raphson inversion and find new increments
    iter = 0
    while (((Jguess - J)^2 + (Lguess - L)^2) > (params.invε)^2)

        Jguess, Lguess, dJgda, dLgda, dJgde, dLgde = ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,d3ψ,aguess,eguess,params)

        # one break: negative actions when getting very close to the centre.
        # define the failure mode: return circular orbit at the minimum size
        if (Jguess < 0.0) | (Lguess < 0.0)
            return da,0.0,iter+1,1.
        end

        # this increment reports occasional failures; why?
        try
            increment1, increment2 = inverse2Dlinear(dJgda,dJgde,dLgda,dLgde,J-Jguess,L-Lguess) 
            aguess, eguess = aguess + increment1, eguess + increment2
        catch e
            # reset to 'safe' values
            aguess, eguess = aguess + 10*da, 0.5
            increment1, increment2 = 0.0, 0.0

            if iter > params.ITERMAX
                finaltol = ((J - Jguess)^2 + (L - Lguess)^2)
                return aguess,eguess,-2,finaltol
            end
        end

        # the try...catch above is failing for some reason
        if (@isdefined increment1) == false
            increment1, increment2 = 0.0, 0.0
        end

        # @WARNING: these appear to have broken something.
        # if bad guesses, needs to reset to a different part of space
        # can't go too small
        if eguess < TOLECC
            # go halfway between the previous guess and 0.
            try
                # reset eguess value
                eguess = eguess - increment2
                eguess = max(TOLECC,0.5eguess)
            catch e
                #(VERBOSE > 1) && println("OrbitalElements.NumericalInversion.AEFromΩ1Ω2Brute: guessing close to ecc=0: ",eguess," (a=",aguess,")")
                eguess = max(TOLECC,0.5eguess)
            end
        end

        if eguess >= (1.0-TOLECC)
            # go halfway between the previous guess and 1.
            try
                eguess = eguess - increment2
                eguess = min(1.0-TOLECC,eguess + 0.5*(1-eguess))
            catch e
                println("OrbitalElements.NumericalInversion.AEFromJLBrute: guessing close to ecc=1: ",eguess," (a=",aguess,") for increment ",increment1," ",increment2)
                eguess = min(1.0-TOLECC,eguess + 0.5*(1-eguess))
            end
        end

        iter += 1

        if iter > params.ITERMAX
            break
        end

    end

    finaltol = ((J - Jguess)^2 + (L - Lguess)^2)

    # check here to not allow bad values?
    if isnan(aguess) | isnan(eguess)
        #(VERBOSE > 1) && println("OrbitalElements.NumericalInversion.AEFromΩ1Ω2Brute: failed for inputs (Ω₁,Ω₂)=($Ω₁,$Ω₂).")

        # return a semi-equivalent circular orbit, as the failure mode is mostly very small orbits
        return acirc,0.0,-1,finaltol

    else
        return aguess,eguess,iter,finaltol
    end
end