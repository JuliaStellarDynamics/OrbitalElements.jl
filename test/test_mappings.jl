# Does not test ...
# - Analytic functions against known values
# - Potentials other than isochrone and Plummer
# @IMPROVE list untested stuff

@testset "mappings" begin
    @testset "model independent" begin
        a, e = 1.0, 0.5
        rp, ra = rpra_from_ae(a, e)
        @test all(rpra_from_ae(1.0, 0.0) .≈ (1.0, 1.0))
        @test rpra_to_ae_jacobian(rp, ra) * ae_to_rpra_jacobian(a, e) ≈ 1
    end
    ###############################
    # Forward mapping test function
    ###############################
    function compare_mappings(
        mapping_analytic::Function,
        mapping_numerical::Function,
        arange,
        erange;
        atol::Float64=0.0,
        rtol::Float64=eps(Float64)
    )
        # Compare the results from analytical and numerical mappings
        for a in arange
            for e in erange
                @test all(
                    isapprox.(
                        mapping_analytic(a, e),
                        mapping_numerical(a, e),
                        atol=atol, 
                        rtol=rtol
                    )
                )
            end
        end
    end
    ################################
    # Backward mapping test function
    ################################
    function test_backwardmapping(
        mapping_forward::Function,
        mapping_backward::Function,
        arange,
        erange,
        tol::Float64
    )
        for a in arange
            for e in erange
                x1, x2 = mapping_forward(a, e)
                @test all(isapprox.(mapping_backward(x1, x2), (a, e), atol=tol, rtol=tol))
            end
        end
    end
    # Range test
    aregular = 0.1:0.1:10.
    aborder = 0.:0.005:0.1
    eregular = 0.1:0.01:0.9
    eborders = Iterators.flatten((0.:0.005:0.1,0.9:0.005:1.0))
    #################################
    # Isochrone Analytic vs Numerical
    #################################
    @testset "Isochrone" begin
        # Compare analytical to numerical results in the isochrone
        anapot = AnalyticIsochrone()
        numpot = NumericalIsochrone()
        anaparams = OrbitalParameters(rc=radial_scale(anapot))
        numparams = OrbitalParameters(rc=radial_scale(numpot))
        @testset "forward" begin
            for (mapping, fun) in [
                ("EL", EL_from_ae),
                ("actions", actions_from_ae),
                ("frequencies", frequencies_from_ae)
            ]
                @testset "$mapping" begin
                    # Defining the mappings to compare
                    ana(a, e) = fun(a, e, anapot, anaparams)
                    num(a, e) = fun(a, e, numpot, numparams)
                    @testset "regular" begin
                        tol =  1.e-6
                        compare_mappings(ana, num, aregular, eregular; atol=tol, rtol=tol)
                    end
                    @testset "borders" begin
                        tol =  1.e-3
                        compare_mappings(ana, num, aborder, eregular; atol=tol, rtol=tol)
                        compare_mappings(ana, num, aregular, eborders; atol=tol, rtol=tol)
                        compare_mappings(ana, num, aborder, eborders; atol=tol, rtol=tol)
                    end
                end
            end
        end
        @testset "backward" begin
            for (version, pot, params) in [
                ("analytic", anapot, anaparams), 
                ("numerical", numpot, numparams)
            ]
                @testset "$version" begin
                    for (mapping, forwardfun, backwardfun) in [
                        ("EL", EL_from_ae, ae_from_EL),
                        ("actions", actions_from_ae, ae_from_actions),
                        ("frequencies", frequencies_from_ae, ae_from_frequencies)
                    ]
                        @testset "$mapping" begin
                            # Defining the mappings to compare
                            forward(a, e) = forwardfun(a, e, pot, params)
                            backward(E, L) = backwardfun(E, L, pot, params)
                            @testset "regular" begin
                                tol =  1.e-6
                                test_backwardmapping(
                                    forward, backward, aregular, eregular, tol
                                )
                            end
                            @testset "borders" begin
                                tol =  1.e-3
                                test_backwardmapping(
                                    forward, backward, aregular, eborders, tol
                                    )
                                # @IMPROVE: for now still issues in the centre !
                                # test_backwardmapping(
                                #     forward, backward, aborder, eregular, tol
                                # )
                                # test_backwardmapping(
                                #     forward, backward, aborder, eborders, tol
                                # )
                            end
                        end
                    end
                end
            end
        end
        @testset "resonant" begin
            # Compare analytical to numerical results in the isochrone
            anapot = AnalyticIsochrone()
            numpot = NumericalIsochrone()
            params = OrbitalParameters(rc=radial_scale(numpot))
            
            @testset "frequency_extrema" begin
                
            end
            numres = Resonance(-1, 2, numpot, numparams)
            # Forward/backward
            a, e = 1.0, 0.5
            α, β = αβ_from_ae(a, e, numpot, numparams)
            # @IMPROVE: to continue !
    
        end
    end
    #################################
    # Plummer Tep vs Numerical
    #################################
    @testset "Plummer" begin
        # Compare analytical to numerical results in the isochrone
        anapot = SemiAnalyticPlummer()
        numpot = NumericalPlummer()
        anaparams = OrbitalParameters(rc=radial_scale(anapot))
        numparams = OrbitalParameters(rc=radial_scale(numpot))
        @testset "forward" begin
            for (mapping, fun) in [
                ("EL", EL_from_ae),
                ("actions", actions_from_ae)
                # ("frequencies", frequencies_from_ae)
                # @IMPROVE: Fail for frequency, for now issues with effective radius
                # and anomaly. Ω_ϕ is wrong when computed with SemiAnalyticPlummer
            ]
                @testset "$mapping" begin
                    # Defining the mappings to compare
                    ana(a, e) = fun(a, e, anapot, anaparams)
                    num(a, e) = fun(a, e, numpot, numparams)
                    @testset "regular" begin
                        tol =  1.e-6
                        compare_mappings(ana, num, aregular, eregular; atol=tol, rtol=tol)
                    end
                    @testset "borders" begin
                        tol =  1.e-3
                        compare_mappings(ana, num, aborder, eregular; atol=tol, rtol=tol)
                        compare_mappings(ana, num, aregular, eborders; atol=tol, rtol=tol)
                        compare_mappings(ana, num, aborder, eborders; atol=tol, rtol=tol)
                    end
                end
            end
        end
        @testset "backward" begin
            for (version, pot, params) in [
                ("analytic", anapot, anaparams), 
                ("numerical", numpot, numparams)
            ]
                @testset "$version" begin
                    for (mapping, forwardfun, backwardfun) in [
                        ("EL", EL_from_ae, ae_from_EL),
                        ("actions", actions_from_ae, ae_from_actions)
                        # ("frequencies", frequencies_from_ae, ae_from_frequencies)
                        # @IMPROVE: Fail for frequency, for now issues with effective radius
                        # and anomaly. Ω_ϕ is wrong when computed with SemiAnalyticPlummer
                    ]
                        @testset "$mapping" begin
                            # Defining the mappings to compare
                            forward(a, e) = forwardfun(a, e, pot, params)
                            backward(E, L) = backwardfun(E, L, pot, params)
                            @testset "regular" begin
                                tol =  1.e-6
                                test_backwardmapping(
                                    forward, backward, aregular, eregular, tol
                                )
                            end
                            @testset "borders" begin
                                tol =  1.e-3
                                test_backwardmapping(
                                    forward, backward, aregular, eborders, tol
                                )
                                # @IMPROVE: for now still issues in the centre !
                                # test_backwardmapping(
                                #     forward, backward, aborder, eregular, tol
                                # )
                                # test_backwardmapping(
                                #     forward, backward, aborder, eborders, tol
                                # )
                            end
                        end
                    end
                end
            end
        end
    end
end