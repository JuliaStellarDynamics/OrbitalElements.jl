"""
test to compare Jacobians

"""


import OrbitalElements

# define the potential to pass to frequency calculators
const bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64   = OrbitalElements.isochrone_psi(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.isochrone_dpsi_dr(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.isochrone_ddpsi_ddr(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.isochrone_dddpsi_dddr(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.isochrone_ddddpsi_ddddr(r,bc,M,G)

Ω₀      = OrbitalElements.isochrone_Omega0(bc,M,G)


# set up an array of semi-major axes to try: log spaced
NORBITS = 1000
x = LinRange(-3.0,-1.0,NORBITS)

function test_freq_derivative_jacobian()
    for o=1:NORBITS
        sma = 10.0^(x[o])
        ecc = 0.4

        # version using numerical derivatives of frequencies
        JacEL1 = OrbitalElements.JacELToAlphaBetaAE(sma,ecc,ψ,dψ,d2ψ,Ω₀,NINT=512)
    end
end

function test_freq_theta_jacobian()
    for o=1:NORBITS
        sma = 10.0^(x[o])
        ecc = 0.4

        # version using theta relationships
        JacEL2 = OrbitalElements.JacELToAlphaBetaAE(ψ,dψ,d2ψ,d3ψ,d4ψ,sma,ecc,Omega0=Ω₀,NINT=64,EDGE=0.03)

    end
end

function test_freq_analytic_jacobian()
    for o=1:NORBITS
        sma = 10.0^(x[o])
        ecc = 0.4

        # the true isochrone value: need alpha, beta first
        Ω₁e,Ω₂e = OrbitalElements.IsochroneOmega12FromAE(sma,ecc,bc,M,G)
        alpha,beta = Ω₁e/Ω₀,Ω₂e/Ω₁e
        JacEL3 = OrbitalElements.IsochroneJacELtoAlphaBeta(alpha,beta,bc,M,G)

        # compute relative errors
        #relerror1 = abs(JacEL1-JacEL3)/JacEL3
        #relerror2 = abs(JacEL2-JacEL3)/JacEL3

        #println("sma=$sma, numderiv=$relerror1, thetaderiv=$relerror2")
        #println("sma=$sma, #1=$JacEL1, #2=$JacEL2, #3=$JacEL3")
    end
end

@time test_freq_derivative_jacobian()
@time test_freq_theta_jacobian()
@time test_freq_analytic_jacobian()
