# Does not test ...
# - Analytic functions against known values
# - Potentials other than isochrone
# @IMPROVE list untested stuff

@testset "mappings" begin
    # Compare analytical to numerical results in the isochrone
    anapot = AnalyticIsochrone()
    numpot = NumericalIsochrone()
    params = OrbitalParameters(Ω₀=Ω₀(numpot))
    tol =  1.e-6
    @testset "forward" begin
        a, e = 1.0, 0.5
        # E, L
        @test all(EL_from_ae(a, e, anapot) .≈ EL_from_ae(a, e, numpot, params))
        # J, L
        @test all(
            isapprox.(
                actions_from_ae(a, e, anapot),
                actions_from_ae(a, e, numpot, params),
                atol=tol
            )
        )
        # Ω₁, Ω₂
        @test all(
            isapprox.(
                frequencies_from_ae(a,e,anapot),
                frequencies_from_ae(a, e, numpot, params),
                atol=tol
            )
        )
    end
    @testset "backward" begin
        a, e = 1.0, 0.5
        tol =  1.e-6
        # E, L
        E, L = EL_from_ae(a, e, anapot)
        @test all(isapprox.(ae_from_EL(E, L, anapot), (a, e), atol=tol))
        @test all(isapprox.(ae_from_EL(E, L, numpot, params), (a, e), atol=tol))
        # J, L
        J, L = actions_from_ae(a, e, anapot)
        @test all(isapprox.(ae_from_actions(J, L, numpot, params), (a, e), atol=tol))
        # Ω₁, Ω₂
        Ω1, Ω2 = frequencies_from_ae(a, e, anapot)
        @test all(isapprox.(ae_from_frequencies(Ω1, Ω2, anapot), (a, e), atol=tol))
        @test all(isapprox.(ae_from_frequencies(Ω1, Ω2, numpot, params), (a, e), atol=tol))
    end
    @testset "resonance" begin 
    end
end