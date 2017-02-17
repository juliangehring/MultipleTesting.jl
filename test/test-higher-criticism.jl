## Test higher criticism ##
module Test_HigherCriticism

using MultipleTesting
using Base.Test


@testset "Higher criticism" begin

    pval1 = [0.00001, 0.0001, 0.002, 0.002, 0.1, 0.1, 0.2]
    hcs1 = [1.080048, 1.672734, 3.044358, 3.044358, 5.724654, 5.724654, 6.048691]
    hcv1 = 0.2

    pval2 = [0, 1e-10, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
    hcs2 = [1.032796, 1.511858, 1.921538, 2.309400, 2.696791, 3.098304, 3.526862, 3.992000,
            3.729249, 3.511505, 3.344031, 3.233162, 3.202563, 3.326087, 2.272150, 0.0]
    hcv2 = 0.001

    @testset "Higher criticism scores" begin

        higher_criticism_scores = MultipleTesting.higher_criticism_scores

        hcs = higher_criticism_scores(PValues(pval1))
        @test isapprox( hcs, hcs1, atol = 1e-6 )

        hcs = higher_criticism_scores(PValues(pval2))
        @test isapprox( hcs, hcs2, atol = 1e-6 )

    end


    @testset "Higher critical value" begin

        higher_critical_value = MultipleTesting.higher_critical_value

        hcv = higher_critical_value(PValues(pval1))
        @test hcv == hcv1

        hcv = higher_critical_value(PValues(pval2))
        @test hcv == hcv2

    end

end

end
