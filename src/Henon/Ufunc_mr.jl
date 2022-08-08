function henon_f(u::Float64)
    return u*(1.5 - 0.5*(u)^(2))
end
function henon_df(u::Float64)
    return 1.5*(1.0 - (u)^(2))
end
function henon_d2f(u::Float64)
    return -3u
end
function henon_d3f(u::Float64)
    return -3.
end
function henon_d4f(u::Float64)
    return 0.
end

function ψeff(ψ::Function,r::Float64,L::Float64)
    return ψ(r) + 0.5 * (L/r)^(2)
end

function dψeffdr(dψ::Function,r::Float64,L::Float64)
    return dψ(r) - (L)^(2) / (r^3)
end

function d2ψeffdr2(d2ψ::Function,r::Float64,L::Float64)
    return d2ψ(r) + 3 * (L)^(2) / (r^4)
end

function ru(fun::Function,u::Float64,a::Float64,e::Float64)
    return a*(1+e*fun(u))
end

function Θ_expansion(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                    u::Float64,
                    a::Float64,e::Float64;
                    TOLECC::Float64=ELTOLECC,
                    f::Function=henon_f,df::Function=henon_df,d2f::Function=henon_d2f,d3f::Function=henon_d3f,d4f::Function=henon_d4f)

    ul = (u > 0) ? 1.0 : -1.0 
    rl = ru(f,ul,a,e)

    E, L = EL_from_ae_pot(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)

    ψeffl, dψeffl, d2ψeffl = ψeff(ψ,rl,L), dψeffdr(dψ,rl,L), d2ψeffdr2(d2ψ,rl,L)

    fl, dfl, d2fl, d3fl, d4fl = f(ul), df(ul), d2f(ul), d3f(ul), d4f(ul)    
    
    pref = - ul * sqrt(2.0) * a * e
    denom = sqrt(- a * e * dψeffl * d2fl)

    a = d2fl
    b = d3fl / 3.0
    c = (
        3.0 * (dψeffl * d2fl * d4fl - a * e * d2ψeffl * (d2fl)^(3)) 
        - dψeffl * (d3fl)^(2)
        ) 
        / (24 * dψeffl * d2fl)

    return pref / den * ( a + b * (u - ul) + c * (u - ul)^(2) )
end


function Θ(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
            u::Float64,
            a::Float64,e::Float64;
            ulim::Float64=0.01,
            TOLECC::Float64=ELTOLECC,
            f::Function=henon_f,df::Function=henon_df,d2f::Function=henon_d2f,d3f::Function=henon_d3f,d4f::Function=henon_d4f)

    if ((1-abs(u))<ulim)     # use the expanded approximation
        return Θ_expansion(ψ,dψ,d2ψ,d3ψ,u,a,e;fun=fun,TOLECC=TOLECC)
    else 
        E, L = EL_from_ae_pot(ψ,dψ,d2ψ,d3ψ,a,e;TOLECC=TOLECC)
        r = ru(f,u,a,e)

        return a * e * df(u) / sqrt(E - ψeff(ψ,r,L)) 
    end
end



function integrate_u(fun,ra::Float64)
    Delta_u = 2.0 / K_MF
    u = -1.0
    res = Delta_u * fun(ra,u) / 6.0

    for k = 1:(K_MF-1)
        u += 0.5*Delta_u
        res += 2.0 * Delta_u * fun(ra,u) / 3.0
        u += 0.5*Delta_u
        res += 1.0 * Delta_u * fun(ra,u) / 3.0
    end
    u += 0.5*Delta_u
    res += 2.0 * Delta_u * fun(ra,u) / 3.0
    u += 0.5*Delta_u
    res += 1.0 * Delta_u * fun(ra,u) / 6.0

    return res
end