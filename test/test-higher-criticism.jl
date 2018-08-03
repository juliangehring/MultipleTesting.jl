## Test higher criticism ##
module Test_HigherCriticism

using MultipleTesting
using Test


@testset "Higher criticism" begin

    pval1 = [0.00001, 0.0001, 0.002, 0.002, 0.1, 0.1, 0.2]
    hcs1_ref = [1.080048, 1.672734, 3.044358, 3.044358, 5.724654, 5.724654, 6.048691]
    hcv1_ref = 0.2

    pval2 = [0, 1e-10, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
    hcs2_ref = [1.032796, 1.511858, 1.921538, 2.309400, 2.696791, 3.098304, 3.526862, 3.992000,
            3.729249, 3.511505, 3.344031, 3.233162, 3.202563, 3.326087, 2.272150, 0.0]
    hcv2_ref = 0.001


    @testset "Higher criticism scores" begin

        hcs1 = estimate(PValues(pval1), HigherCriticismScores())
        @test isapprox( hcs1, hcs1_ref, atol = 1e-5 )

        hcs2 = estimate(PValues(pval2), HigherCriticismScores())
        @test isapprox( hcs2, hcs2_ref, atol = 1e-5 )

        @test_throws MethodError estimate(PValues(pval1))
        @test_throws MethodError estimate(HigherCriticismScores())

        # currently not supported
        @test_throws MethodError estimate(pval1, HigherCriticismScores())
        @test_throws MethodError estimate(PValues(pval1), HigherCriticismScores)

    end


    @testset "Higher criticism threshold" begin

        hcv1 = estimate(PValues(pval1), HigherCriticismThreshold())
        @test hcv1 == hcv1_ref

        hcv2 = estimate(PValues(pval2), HigherCriticismThreshold())
        @test hcv2 == hcv2_ref

        @test_throws MethodError estimate(PValues(pval1))
        @test_throws MethodError estimate(HigherCriticismThreshold())

        # currently not supported
        @test_throws MethodError estimate(pval1, HigherCriticismThreshold())
        @test_throws MethodError estimate(PValues(pval1), HigherCriticismThreshold)

    end

end

end
