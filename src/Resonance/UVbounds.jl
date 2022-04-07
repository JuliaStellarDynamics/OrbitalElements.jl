

function find_w_min_max(n1::Int64,n2::Int64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.,Omega0::Float64=1.)


    extreme(x) = n1*Omega1_circular(dpotential,ddpotential,x) + n2*Omega2_circular(dpotential,x)

    m = extremise_function(extreme,24,0.,rmax,false) # return a min/max flag?

    # for the problem of a cored cluster with an infinite extent, the w minima and maxima are either
    # in the very centre
    # in the very outskirts
    # along the circular velocity
    w_min,w_max = extrema([extreme(1.e-10)/Omega0,0.0,extreme(m)/Omega0])
    return w_min,w_max
end


function get_varpi(omg::Complex{Float64},n1::Int64,n2::Int64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.,Omega0::Float64=1.)


    w_min,w_max = find_w_min_max(n1,n2,dpotential,ddpotential,rmax,Omega0)

    return (2omg - w_max - w_min)/(w_max - w_min)

end


#w1,w2 = find_w_min_max(3,-4,isochrone_dpsi_dr,isochrone_ddpsi_ddr)
#print(w1," ",w2,"\n")



function make_betac(dpotential::Function,ddpotential::Function,numr::Int64=2000,Omega0::Float64=1.)
    # in epicycle.jl
    # do a high-resolution interpolation to get \beta_c(alpha)

    alpha_c(x) = Omega1_circular(dpotential,ddpotential,x)    # alpha_c(r)
    beta_c(x)  = Omega2_circular(dpotential,x)/alpha_c(x)                   # beta_c(r)
    # so we need to invert g(x) to find x, then use to solve f (when n2!=0)
    # when n2 = 0, we need beta_c as a function of beta

    testu = 10 .^ LinRange(5.,-5.,numr)

    garr = Array{Float64}(undef, (numr))
    farr = Array{Float64}(undef, (numr))

    for u = 1:numr
        garr[u] = alpha_c(testu[u])/Omega0
        farr[u] = beta_c(testu[u])
    end

    beta_c = LinearInterpolation(garr,farr)

    return beta_c
end


function hu(u::Float64,wmin::Float64,wmax::Float64)
    return 0.5*(wmax+wmin + u*(wmax-wmin))
end


function root_of_h_omega(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64,vbound::Float64,beta_c)
    #=
    requires beta_c, should pass this.
    =#
    hval = hu(u,wmin,wmax)
    rootequation(x) = hval - n1*x - n2*x*beta_c[x]
    # two roots to try: bounded by [0,v(u=1)] and [v(u=1),1]
    if rootequation(0+1.e-6)*rootequation(vbound) < 0
        r1 = fzero(rootequation, 0+1.e-6,vbound)
    else
        r1 = 0.0
    end

    if rootequation(1-1.e-6)*rootequation(vbound) < 0
        r2 = fzero(rootequation, vbound,1-1.e-6)
    else
        r2 = 1.0
    end
    #print(r1," ",r2)

    # edge curing for vbounds
    # from upper left to lower right cure
    if (vbound > 0.999999) & (u==-1.) & (hu(-1.,wmin,wmax)<0)
        r1 = 0.999999
        r2 = 1.0
    end

    # from lower left to upper right cure
    if (vbound > 0.999999) & (u==1.) & (hu(1.,wmin,wmax)>0)
        r1 = 0.999999
        r2 = 1.
    end

    # greater than r1, less than r2
    return r1,r2
end

function constraint_three(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64)
    # greater than this
    hval = hu(u,wmin,wmax)
    return hval/(n2/2 + n1)
end

function vmin_vmax(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64,vbound::Float64,beta_c)
    # this function works for n2 != 0

    if (n2==0)

        hval = hu(u,wmin,wmax)

        vmin = 0.5

        # put in guards for the very edges. SLOPPY
        vmax = beta_c(minimum([0.99999,maximum([hval/n1,0.00001])]))

    else

        r1,r2 = root_of_h_omega(u,wmin,wmax,n1,n2,vbound,beta_c)

        r3 = constraint_three(u,wmin,wmax,n1,n2)

        if (abs(n2)<=abs(n1)) | (abs(n2)>=abs(2n1))
            vmin = maximum([0.0,r1])
            vmax = minimum([r2,1.0,r3])
        else
            vmin = maximum([0.0,r1,r3])
            vmax = minimum([r2,1.0])
        end

        if isnan(vmax)
            vmax = 1.
        end

        if isnan(vmin)
            vmin = 0.
        end
    end

    return vmin,vmax
end

# test this transformation, set up the bilinear interpolation in alpha, beta

function alphabeta_from_uv(u::Float64,v::Float64,
                           n1::Int64,n2::Int64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.)
    #=
    the inverse mapping from (u,v) -> (alpha,beta)

    weirdly imperfect: why?

    =#

    wmin,wmax = find_w_min_max(n1,n2,dpotential,ddpotential,rmax)

    if n2 == 0
        beta  = v
        alpha = (1/(2n1))*((wmax-wmin)*u + wmin + wmax)
    else
        alpha = v
        beta  = (1/(n2*v))*(0.5*((wmax-wmin)*u + wmin + wmax) - n1*v)
    end

    return alpha,beta
end

function uv_from_alphabeta(alpha::Float64,beta::Float64,
                           n1::Int64,n2::Int64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.)
    #=
    the mapping from (alpha,beta) -> (u,v)

    =#
    wmin,wmax = find_w_min_max(n1,n2,dpotential,ddpotential,rmax)

    wval = n1*alpha + n2*beta*alpha

    u = (2*wval - wmax - wmin)/(wmax-wmin)

    if (n2==0)
        v = beta
    else
        v = alpha
    end

    return u,v

end
