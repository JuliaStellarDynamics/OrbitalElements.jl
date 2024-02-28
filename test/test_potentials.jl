# Does not test ...
# - derivatives (against known values or finite differences ?)
# - Scale functions (E₀,L₀,Ω₀)
# @IMPROVE list untested stuff

@testset "potentials" begin
    # # Check functions
    # function checkapprox(f,x,y;rtol::Real=√eps)
    #     # Verify input/output compatibility
    #     @assert length(x) == length(y)
    #     nchecks = length(x) # Number of checked values
    #     # Checks
    #     for i=1:nchecks; @test f(x[i]) ≈ y[i] rtol=rtol; end
    # end
    # function checkdomainerror(f,x)
    #     for locx in x; @test_throws DomainError f(locx); end
    # end
    # # Check locations
    # rcheck = (0., 1., Inf) # Radii to check for values
    # rerr = (-1.) # Radii with expected (domain) error
    @testset "isochrone" begin
        # Isochrone potential with default values
        # supposedly G=1, M=1, bc=1
        pot = IsochronePotential()
        # Check potential values
        @test_throws DomainError ψ(pot,-1.)
        @test ψ(pot,0.) ≈ -0.5
        @test ψ(pot,1.) ≈ 1-sqrt(2)
        @test ψ(pot,Inf) == 0.
        # Check derivatives
        # Check scale functions
    end
    @testset "plummer" begin
        # Plummer potential with default values
        # supposedly G=1, M=1, bc=1
        pot = PlummerPotential()
        @test_throws DomainError ψ(pot,-1.)
        @test ψ(pot,0.) ≈ -1.
        @test ψ(pot,1.) ≈ -1/sqrt(2)
        @test ψ(pot,Inf) == 0.
        # Check derivatives
        # Check scale functions
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
            # Check derivatives
            # Check scale functions
        end
        @testset "tapered" begin
            # Tapered Mestel potential with default values
            # supposedly R0=1, V0=1, ε0=1.e-5
            pot = TaperedMestel()
            @test_throws DomainError ψ(pot,-1.)
            @test ψ(pot,0.) ≈ -5log(10)
            @test ψ(pot,1.) ≈ 5.e-11 rtol=1.e-5
            @test ψ(pot,Inf) == Inf
            # Check derivatives
            # Check scale functions
        end
    end
end