# Does not test ...
# @IMPROVE list untested stuff

# @IMPORTANT we rely on potential being defined inside this @testset begin ... end (i.e. inside a "run" function) to avoid allocations
@testset "allocs" begin
    anapot = AnalyticIsochrone()
    numpot = NumericalIsochrone()
    params = OrbitalParameters(rc=radial_scale(numpot))
    tested_ae = [(0.001, 0.001), (0.5, 0.999), (10., 0.001), (1.0, 0.5)]
    @testset "forward" begin
        for (a, e) in tested_ae
            # E, L
            @test (@allocated EL_from_ae(a, e, anapot)) == 0
            @test (@allocated EL_from_ae(a, e, anapot, params)) == 0
            @test (@allocated EL_from_ae(a, e, numpot)) == 0
            @test (@allocated EL_from_ae(a, e, numpot, params)) == 0
            # J, L
            @test (@allocated actions_from_ae(a, e, anapot)) == 0
            @test (@allocated actions_from_ae(a, e, anapot, params)) == 0
            @test (@allocated actions_from_ae(a, e, numpot)) == 0
            @test (@allocated actions_from_ae(a, e, numpot, params)) == 0
            # Ω₁, Ω₂
            @test (@allocated frequencies_from_ae(a, e, anapot)) == 0
            @test (@allocated frequencies_from_ae(a, e, anapot, params)) == 0
            @test (@allocated frequencies_from_ae(a, e, numpot)) == 0
            @test (@allocated frequencies_from_ae(a, e, numpot, params)) == 0
        end
    end
    @testset "backward" begin
        for (a, e) in tested_ae
            # E, L
            E, L = EL_from_ae(a, e, anapot)
            @test (@allocated ae_from_EL(E, L, anapot)) == 0
            @test (@allocated ae_from_EL(E, L, anapot, params)) == 0
            @test (@allocated ae_from_EL(E, L, numpot)) == 0
            @test (@allocated ae_from_EL(E, L, numpot, params)) == 0
            # J, L
            J, L = actions_from_ae(a, e, anapot)
            @test (@allocated ae_from_actions(J, L, anapot)) == 0
            @test (@allocated ae_from_actions(J, L, anapot, params)) == 0
            @test (@allocated ae_from_actions(J, L, numpot)) == 0
            @test (@allocated ae_from_actions(J, L, numpot, params)) == 0
            # Ω₁, Ω₂
            Ω₁, Ω₂ = frequencies_from_ae(a, e, anapot)
            @test (@allocated ae_from_frequencies(Ω₁, Ω₂, anapot)) == 0
            @test (@allocated ae_from_frequencies(Ω₁, Ω₂, anapot, params)) == 0
            @test (@allocated ae_from_frequencies(Ω₁, Ω₂, numpot)) == 0
            @test (@allocated ae_from_frequencies(Ω₁, Ω₂, numpot, params)) == 0
        end
    end
    @testset "resonance" begin 
    end
end