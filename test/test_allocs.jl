# Does not test ...
# @IMPROVE list untested stuff

# @IMPORTANT we rely on potential being defined inside this @testset begin ... end (i.e. inside a "run" function) to avoid allocations
@testset "allocs" begin
    pot = IsochronePotential()
    params = OrbitalParameters(Ω₀=Ω₀(pot))
    @testset "forward" begin
        a, e = 1.0, 0.5
        # E, L
        @test (@allocated ELFromAE(pot,a,e)) == 0
        @test (@allocated ELFromAE(pot,a,e,params)) == 0
        # J, L
        @test (@allocated ComputeActionsAE(pot,a,e)) == 0
        @test (@allocated ComputeActionsAE(pot,a,e,params)) == 0
        # Ω₁, Ω₂
        @test (@allocated ComputeFrequenciesAE(pot,a,e)) == 0
        @test (@allocated ComputeFrequenciesAE(pot,a,e,params)) == 0
    end
    @testset "backward" begin
        a, e = 1.0, 0.5
        # E, L
        E, L = ELFromAE(pot,a,e)
        @test (@allocated AEFromEL(pot,E,L)) == 0
        @test (@allocated OrbitalElements.AEFromELBrute(E,L,pot,params)) == 0
        # J, L
        J, L = ComputeActionsAE(pot,a,e)
        @test (@allocated ComputeAEFromActions(pot,a,e)) == 0
        @test (@allocated ComputeAEFromActions(pot,a,e,params)) == 0
        # Ω₁, Ω₂
        Ω₁, Ω₂ = ComputeFrequenciesAE(pot,a,e)
        @test (@allocated ComputeAEFromFrequencies(pot,Ω₁,Ω₂)) == 0
        @test (@allocated ComputeAEFromFrequencies(pot,Ω₁,Ω₂,params)) == 0
    end
    @testset "resonance" begin 
    end
end