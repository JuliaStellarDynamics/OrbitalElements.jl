using SpecialFunctions

#####
#
#   The Kalnajs' isochrone disc DF
#
#####

"""
    IsoDiscKal_DF(E, L)
Kalnajs distribution function for isochrone disc.
"""
function IsoDiscKal_DF(E::Float64,L::Float64)

    return 0.0
end
"""
    IsoDiscKal_dDFdE(E, L)
Kalnajs DF derivative w.r.t. E.
"""
function IsoDiscKal_dDFdE(E::Float64,L::Float64)

    return 0.0
end
"""
    IsoDiscKal_dDFdL(E, L)
Kalnajs DF derivative w.r.t. E.
"""
function IsoDiscKal_dDFdL(E::Float64,L::Float64)

    return 0.0
end

"""
    IsoDiscKal_ndDFdJ(n1, n2, E, L, ndotOmega)
Kalnajs DF derivative w.r.t. the actions J.
"""
function IsoDiscKal_ndDFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64)

    dDFdE = IsoDiscKal_dDFdE(E,L)
    dDFdL = IsoDiscKal_dDFdL(E,L)
    
    return ndotOmega*dDFdE + n2*dDFdL
end

#####
#
#   The Pichon & Lynden-Bell's isochrone disc DF
#
#####

"""
    IsoDiscPLB_DF(E, L)
Kalnajs distribution function for isochrone disc.
"""
function IsoDiscPLB_DF(E::Float64,L::Float64)

    return 0.0
end
"""
    IsoDiscPLB_dDFdE(E, L)
Kalnajs DF derivative w.r.t. E.
"""
function IsoDiscPLB_dDFdE(E::Float64,L::Float64)

    return 0.0
end
"""
    IsoDiscPLB_dDFdL(E, L)
Kalnajs DF derivative w.r.t. E.
"""
function IsoDiscPLB_dDFdL(E::Float64,L::Float64)

    return 0.0
end

"""
    IsoDiscPLB_ndDFdJ(n1, n2, E, L, ndotOmega)
Kalnajs DF derivative w.r.t. the actions J.
"""
function IsoDiscPLB_ndDFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64)

    dDFdE = IsoDiscPLB_dDFdE(E,L)
    dDFdL = IsoDiscPLB_dDFdL(E,L)
    
    return ndotOmega*dDFdE + n2*dDFdL
end