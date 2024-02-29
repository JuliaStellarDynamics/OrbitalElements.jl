# Does not test :
#   - undefined scale functions (E₀, L₀ for Mestel discs)
# @IMPROVE add undefined scale functions tests once defined

@testset "potentials" begin
    @testset "isochrone" begin
        # Wrong model characteristic values
        @test_throws DomainError IsochronePotential(M=-1.)
        @test_throws DomainError IsochronePotential(bc=0.)
        # Isochrone potential with default values
        # supposedly G=1, M=1, bc=1
        pot = IsochronePotential()
        # Check potential values
        @test_throws DomainError ψ(pot,-1.)
        @test ψ(pot,0.) ≈ -0.5
        @test ψ(pot,1.) ≈ 1-sqrt(2)
        @test ψ(pot,Inf) == 0.
        # Check first derivative
        @test_throws DomainError dψ(pot,-1.)
        @test dψ(pot,0.) == 0.
        @test dψ(pot,1.) ≈ 3/sqrt(2)-2
        @test dψ(pot,Inf) == 0.
        # Check second derivative
        @test_throws DomainError d2ψ(pot,-1.)
        @test d2ψ(pot,0.) == 0.25
        @test d2ψ(pot,1.) ≈ 6-17sqrt(2)/4
        @test d2ψ(pot,Inf) == 0.
        # Check scale functions
        @test Ω₀(pot) ≈ 1.
        @test E₀(pot) ≈ -1.
        @test L₀(pot) ≈ 1.
    end
    @testset "plummer" begin
        # Wrong model characteristic values
        @test_throws DomainError PlummerPotential(M=-1.)
        @test_throws DomainError PlummerPotential(bc=0.)
        # Plummer potential with default values
        # supposedly G=1, M=1, bc=1
        pot = PlummerPotential()
        @test_throws DomainError ψ(pot,-1.)
        @test ψ(pot,0.) ≈ -1.
        @test ψ(pot,1.) ≈ -1/sqrt(2)
        @test ψ(pot,Inf) == 0.
        # Check first derivative
        @test_throws DomainError dψ(pot,-1.)
        @test dψ(pot,0.) == 0.
        @test dψ(pot,1.) ≈ sqrt(2)/4
        @test dψ(pot,Inf) == 0.
        # Check second derivative
        @test_throws DomainError d2ψ(pot,-1.)
        @test d2ψ(pot,0.) == 1.
        @test d2ψ(pot,1.) ≈ -sqrt(2)/8
        @test d2ψ(pot,Inf) == 0.
        # Check scale functions
        @test Ω₀(pot) ≈ 2.
        @test E₀(pot) ≈ -1.
        @test L₀(pot) ≈ 1.
    end
    @testset "mestel" begin
        @testset "true" begin
            # Mestel potential with default values
            # supposedly R0=1, V0=1
            pot = MestelPotential()
            @test_throws DomainError ψ(pot,-1.)
            @test ψ(pot,0.) == -Inf
            @test ψ(pot,1.) == 0.
            @test ψ(pot,Inf) == Inf
            # Check first derivative
            @test_throws DomainError dψ(pot,-1.)
            @test dψ(pot,0.) == Inf
            @test dψ(pot,1.) ≈ 1.
            @test dψ(pot,Inf) == 0
            # Check second derivative
            @test_throws DomainError d2ψ(pot,-1.)
            @test d2ψ(pot,0.) == -Inf
            @test d2ψ(pot,1.) ≈ -1.
            @test d2ψ(pot,Inf) == 0.
            # Check scale functions
            @test Ω₀(pot) ≈ 1.
        end
        @testset "tapered" begin
            # Tapered Mestel potential with default values
            # supposedly R0=1, V0=1, ε0=1.e-5
            pot = TaperedMestel()
            @test_throws DomainError ψ(pot,-1.)
            @test ψ(pot,0.) ≈ -5log(10)
            @test ψ(pot,1.) ≈ 5.e-11 rtol=1.e-5
            @test ψ(pot,Inf) == Inf
            # Check first derivative
            @test_throws DomainError dψ(pot,-1.)
            @test dψ(pot,0.) == 0.
            @test dψ(pot,1.) ≈ 1. rtol=1.e-5
            @test dψ(pot,Inf) == 0.
            # Check second derivative
            @test_throws DomainError d2ψ(pot,-1.)
            @test d2ψ(pot,0.) ≈ 1.e10 rtol=1.e-5
            @test d2ψ(pot,1.) ≈ -1.
            @test d2ψ(pot,Inf) == 0.
            # Check scale functions
            @test Ω₀(pot) ≈ 1.
        end
    end
end