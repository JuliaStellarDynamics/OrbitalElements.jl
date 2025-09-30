# Does not test :
#   - undefined scale functions (E₀, L₀ for Mestel discs)
# @IMPROVE add undefined scale functions tests once defined

@testset "potentials" begin
    @testset "isochrone" begin
        @testset "analytical" begin
            # Wrong model characteristic values
            @test_throws DomainError AnalyticIsochrone(M=-1.)
            @test_throws DomainError AnalyticIsochrone(bc=0.)
            # Isochrone potential with default values
            # supposedly G=1, M=1, bc=1
            model = AnalyticIsochrone()
            # Check potential values
            @test_throws DomainError ψ(-1.,model)
            @test ψ(0.,model) ≈ -0.5
            @test ψ(1.,model) ≈ 1-sqrt(2)
            @test ψ(Inf,model) == 0.
            # Check first derivative
            @test_throws DomainError dψ(-1.,model)
            @test dψ(0.,model) == 0.
            @test dψ(1.,model) ≈ 3/sqrt(2)-2
            @test dψ(Inf,model) == 0.
            # Check second derivative
            @test_throws DomainError d2ψ(-1.,model)
            @test d2ψ(0.,model) == 0.25
            @test d2ψ(1.,model) ≈ 6-17sqrt(2)/4
            @test d2ψ(Inf,model) == 0.
            # Check scale functions
            @test frequency_scale(model) ≈ 1.
            @test energy_scale(model) ≈ -1.
            @test momentum_scale(model) ≈ 1.
            @test radial_scale(model) ≈ 1.
        end
        @testset "numerical" begin
            # Wrong model characteristic values
            @test_throws DomainError NumericalIsochrone(M=-1.)
            @test_throws DomainError NumericalIsochrone(bc=0.)
            # Isochrone potential
            @test_nowarn NumericalIsochrone()
            # Do not test the potential and derivatives as they go throught the exact 
            # same functions as the analytical version.
        end
    end
    @testset "plummer" begin
        @testset "numerical" begin
            # Wrong model characteristic values
            @test_throws DomainError NumericalPlummer(M=-1.)
            @test_throws DomainError NumericalPlummer(bc=0.)
            # Plummer potential with default values
            # supposedly G=1, M=1, bc=1
            model = NumericalPlummer()
            @test_throws DomainError ψ(-1.,model)
            @test ψ(0.,model) ≈ -1.
            @test ψ(1.,model) ≈ -1/sqrt(2)
            @test ψ(Inf,model) == 0.
            # Check first derivative
            @test_throws DomainError dψ(-1.,model)
            @test dψ(0.,model) == 0.
            @test dψ(1.,model) ≈ sqrt(2)/4
            @test dψ(Inf,model) == 0.
            # Check second derivative
            @test_throws DomainError d2ψ(-1.,model)
            @test d2ψ(0.,model) == 1.
            @test d2ψ(1.,model) ≈ -sqrt(2)/8
            @test d2ψ(Inf,model) == 0.
            # Check scale functions
            @test frequency_scale(model) ≈ 2.
            @test energy_scale(model) ≈ -1.
            @test momentum_scale(model) ≈ 1.
            @test radial_scale(model) ≈ 1.
        end
        @testset "semi-analytic" begin
            # Wrong model characteristic values
            @test_throws DomainError SemiAnalyticPlummer(M=-1.)
            @test_throws DomainError SemiAnalyticPlummer(bc=0.)
            # Isochrone potential
            @test_nowarn SemiAnalyticPlummer()
            # Do not test the potential and derivatives as they go throught the exact 
            # same functions as the analytical version.
        end
    end
    @testset "mestel" begin
        @testset "true" begin
            # Mestel potential with default values
            # supposedly R0=1, V0=1
            model = MestelPotential()
            @test_throws DomainError ψ(-1.,model)
            @test ψ(0.,model) == -Inf
            @test ψ(1.,model) == 0.
            @test ψ(Inf,model) == Inf
            # Check first derivative
            @test_throws DomainError dψ(-1.,model)
            @test dψ(0.,model) == Inf
            @test dψ(1.,model) ≈ 1.
            @test dψ(Inf,model) == 0
            # Check second derivative
            @test_throws DomainError d2ψ(-1.,model)
            @test d2ψ(0.,model) == -Inf
            @test d2ψ(1.,model) ≈ -1.
            @test d2ψ(Inf,model) == 0.
            # Check scale functions
            @test frequency_scale(model) ≈ 1.
            @test radial_scale(model) ≈ 1.
        end
        @testset "tapered" begin
            # Tapered Mestel potential with default values
            # supposedly R0=1, V0=1, ε0=1.e-5
            model = TaperedMestel()
            @test_throws DomainError ψ(-1.,model)
            @test ψ(0.,model) ≈ -5log(10)
            @test ψ(1.,model) ≈ 5.e-11 rtol=1.e-5
            @test ψ(Inf,model) == Inf
            # Check first derivative
            @test_throws DomainError dψ(-1.,model)
            @test dψ(0.,model) == 0.
            @test dψ(1.,model) ≈ 1. rtol=1.e-5
            @test dψ(Inf,model) == 0.
            # Check second derivative
            @test_throws DomainError d2ψ(-1.,model)
            @test d2ψ(0.,model) ≈ 1.e10 rtol=1.e-5
            @test d2ψ(1.,model) ≈ -1.
            @test d2ψ(Inf,model) == 0.
            # Check scale functions
            @test frequency_scale(model) ≈ 1.
            @test radial_scale(model) ≈ 1.
        end
    end
    @testset "hernquist" begin
        @testset "numerical" begin
            # Wrong model characteristic values
            @test_throws DomainError NumericalHernquist(M=-1.)
            @test_throws DomainError NumericalHernquist(bc=0.)
            # Plummer potential with default values
            # supposedly G=1, M=1, bc=1
            model = NumericalHernquist()
            @test_throws DomainError ψ(-1.,model)
            @test ψ(0.,model) ≈ -1.
            @test ψ(1.,model) ≈ -0.5
            @test ψ(Inf,model) == 0.
            # Check first derivative
            @test_throws DomainError dψ(-1.,model)
            @test dψ(0.,model) == 1.0
            @test dψ(1.,model) ≈ 0.25
            @test dψ(Inf,model) == 0.
            # Check second derivative
            @test_throws DomainError d2ψ(-1.,model)
            @test d2ψ(0.,model) == -2.
            @test d2ψ(1.,model) ≈ -0.25
            @test d2ψ(Inf,model) == 0.
            # Check scale functions
            @test frequency_scale(model) ≈ 1.
            @test energy_scale(model) ≈ -1.
            @test momentum_scale(model) ≈ 1.
            @test radial_scale(model) ≈ 1.
        end
    end
end