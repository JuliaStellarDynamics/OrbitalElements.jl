using HypergeometricFunctions

#####
#
#   The Miyamoto DFs for Kuzmin-Toomre disc 
#
#####

"""
    Miyamoto_DF(E, L)
Miyamoto distribution function for Kuzmin-Toomre disc.
"""
function Miyamoto_DF(E::Float64,L::Float64;mM::Int64=1)

    return (2*mM+3) / (2*(pi^2)) * (E^(2*mM+2)) * _₂F₁(-mM, -2-2*mM, 1/2, L^2/(2*E))
end
"""
    Miyamoto_dDFdE(E, L)
Miyamoto DF derivative w.r.t. E.
"""
function Miyamoto_dDFdE(E::Float64,L::Float64;mM::Int64=1)

    return ((mM+1)*(2*mM+3) / (pi^2)) * (E^(2*mM)) * 
                ( E * _₂F₁(-mM,-2-2*mM,1/2,L^2/(2*E))
                - (L^2) * _₂F₁(-1-2*mM,1-mM,3/2,L^2/(2*E)) )
end
"""
    Miyamoto_dDFdL(E, L)
Miyamoto DF derivative w.r.t. L.
"""
function Miyamoto_dDFdL(E::Float64,L::Float64;mM::Int64=1)

    return (2*mM*(mM+1)*(2*mM+3) / (pi^2)) * L * (E^(2*mM+1)) * _₂F₁(-1-2*mM,1-mM,3/2,L^2/(2*E))
end

"""
    Miyamoto_ndDFdJ(n1, n2, E, L, ndotOmega)
Miyamoto DF derivative w.r.t. the actions J.
"""
function Miyamoto_ndDFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64;mM::Int64=1)

    dDFdE = Miyamoto_dDFdE(E,L;mM=mM)
    dDFdL = Miyamoto_dDFdL(E,L;mM=mM)
    
    return ndotOmega*dDFdE + n2*dDFdL
end
