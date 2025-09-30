
"""
    u_v_from_R_z(R::Float64, z::Float64, Delta::Float64) -> Tuple{Float64, Float64}

Calculate the radial and vertical coordinates (u,v) from the spheroidal coordinates (R,z).

See BT08 eq. 3.242 (pg. 227)

# Arguments
- `R::Float64`: The cylindrical radius.
- `z::Float64`: The vertical height.
- `Delta::Float64`: The spheroidal scale parameter.

# Returns
- `Tuple{Float64, Float64}`: A tuple containing the radial coordinate `u` and the vertical coordinate `v`.
"""
function u_v_from_R_z(R::Float64, z::Float64, Delta::Float64)

    zDeltaplus = z+Delta
    zDeltaminus = z-Delta
    term1 = sqrt(R^2+zDeltaplus^2)
    term2 = sqrt(R^2+zDeltaminus^2)
    prefac = 1/(2Delta)

    u = acosh(prefac*(term1+term2))
    v = acos(prefac*(term1-term2))


    return u,v
end 

"""
    R_z_from_u_v(u::Float64, v::Float64, Delta::Float64) -> Tuple{Float64, Float64}

Calculate the radial and vertical coordinates (R, z) from the parameters (u, v).

# Arguments
- `u::Float64`: The u parameter.
- `v::Float64`: The v parameter.

# Returns
- `Tuple{Float64, Float64}`: A tuple containing the radial coordinate `R` and the vertical coordinate `z`.
"""
function R_z_from_u_v(u::Float64, v::Float64, Delta::Float64)
    R = Delta*sinh(u)*sin(v)
    z = Delta*cosh(u)*cos(v)

    return R, z
end


"""
    R_z_from_lambda_nu(lambda::Float64, nu::Float64) -> Tuple{Float64, Float64}

Calculate the radial and vertical coordinates (R, z) from the spheroidal coordinates (lambda, nu).

# Arguments
- `lambda::Float64`: The lambda coordinate.
- `nu::Float64`: The nu coordinate.

# Returns
- `Tuple{Float64, Float64}`: A tuple containing the radial coordinate `R` and the vertical coordinate `z`.
"""
function R_z_from_lambda_nu(lambda::Float64, nu::Float64, alpha::Float64, gamma::Float64)
    RSq = (lambda+alpha)*(nu+alpha)/(alpha-gamma)
    zSq = (lambda+gamma)*(nu+gamma)/(gamma-alpha)

    return sqrt(RSq), sqrt(zSq)
end 

"""
    lambda_nu_from_u_v(u::Float64, v::Float64) -> Tuple{Float64, Float64}

Calculate the spheroidal coordinates (lambda, nu) from the parameters (u, v).

# Arguments
- `u::Float64`: The u parameter.
- `v::Float64`: The v parameter.

# Returns
- `Tuple{Float64, Float64}`: A tuple containing the lambda coordinate `lambda` and the nu coordinate `nu`.
"""
function lambda_nu_from_u_v(u::Float64, v::Float64, a::Float64, c::Float64)
    lambda = a^2*cosh(u)^2 - c^2*sinh(u)^2
    nu = c^2*sin(v)^2 + a^2*cos(v)^2

    return lambda, nu 
end 

"""
    lambda_nu_from_R_z(u::Float64, v::Float64) -> Tuple{Float64, Float64}

Calculate the spheroidal coordinates (lambda, nu) from the parameters (u, v).

# Arguments
- `R::Float64`: The R parameter.
- `z::Float64`: The z parameter.

# Returns
- `Tuple{Float64, Float64}`: A tuple containing the lambda coordinate `lambda` and the nu coordinate `nu`.
"""
function lambda_nu_from_R_z(R::Float64, z::Float64, a::Float64, c::Float64)

    Delta = sqrt(a^2-c^2)
    u,v = u_v_from_R_z(R,z,Delta)


    lambda, nu = lambda_nu_from_u_v(u,v,a,c)

    return lambda, nu 
end 


"""
    lambda_from_u(u::Float64) -> Float64

Calculate the lambda coordinate from the parameter u.

# Arguments
- `u::Float64`: The u parameter.

# Returns
- `Float64`: The lambda coordinate.
"""
function lambda_from_u(u::Float64, a::Float64, Delta::Float64)

    # lambda = a^2*cosh(u)^2 - c^2*sinh(u)^2
    lambda = a^2+ Delta^2*sinh(u)^2


    return lambda
end 

"""
    nu_from_v(v::Float64) -> Float64

Calculate the nu coordinate from the parameter v.

# Arguments
- `v::Float64`: The v parameter.

# Returns
- `Float64`: The nu coordinate.
"""
function nu_from_v(v::Float64, a::Float64, c::Float64)


    nu = c^2*sin(v)^2 + a^2*cos(v)^2

    return nu 
end 



"""
    xi_eta_from_u_v(u::Float64, v::Float64) -> Tuple{Float64, Float64}

Calculate the xi and eta coordinates from the parameters (u, v).

# Arguments
- `u::Float64`: The u parameter.
- `v::Float64`: The v parameter.

# Returns
- `Tuple{Float64, Float64}`: A tuple containing the xi coordinate `xi` and the eta coordinate `eta`.
"""
function xi_eta_from_u_v(u::Float64, v::Float64)

    xi = cosh(u)
    eta = cos(v)

    return xi, eta
end





