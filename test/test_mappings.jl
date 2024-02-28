# Does not test ...
# - Analytic functions against known values
# - Potentials other than isochrone
# @IMPROVE list untested stuff

@testset "mappings" begin
    # Compare analytical to numerical results in the isochrone
    pot = IsochronePotential()
    params = OrbitalParameters(Ω₀=Ω₀(pot))
    tol =  1.e-6
    @testset "forward" begin
        a, e = 1.0, 0.5
        # E, L
        @test EFromAE(pot,a,e) ≈ EFromAE(pot,a,e,params)
        @test LFromAE(pot,a,e) ≈ LFromAE(pot,a,e,params)
        @test all(ELFromAE(pot,a,e) .≈ ELFromAE(pot,a,e,params))
        # J, L
        @test all(isapprox.(ComputeActionsAE(pot,a,e),ComputeActionsAE(pot,a,e,params),atol=tol))
        # Ω₁, Ω₂
        @test all(isapprox.(ComputeFrequenciesAE(pot,a,e),ComputeFrequenciesAE(pot,a,e,params),atol=tol))
    end
    @testset "backward" begin 
    end
    @testset "resonance" begin 
    end
end