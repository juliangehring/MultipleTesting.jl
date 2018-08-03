## Test pvalue adjustment methods ##
module Test_pval_pi0_adjustment

using MultipleTesting
using Test


@testset "p-Value π₀ adjustment" begin

    pval1 = [0.0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.75, 1.0]
    pi0 = 0.4

    @test_throws MethodError adjust()
    @test_throws MethodError adjust(pval1)


    @testset "benjamini_hochberg" begin

        t = BenjaminiHochbergAdaptive

        ref0 = [0.0, 0.0005, 0.003333333, 0.025, 0.1, 0.166666667, 0.285714286, 0.5, 0.833333333, 1.0]
        ref = ref0 .* pi0

        ## no integers as input
        @test_throws MethodError adjust([0, 1], t())

        ## no valid p-values as input
        @test_throws DomainError adjust([-1.0, 0.5], t(pi0))
        @test_throws DomainError adjust([0.5, 1.5], t(pi0))
        @test_throws DomainError adjust([0.5, 0.7], t(-1.0))
        @test_throws DomainError adjust([0.5, 0.7], t(1.5))

        ## single p-value is returned unchanged
        pval = rand(1)
        @test adjust(pval, t(pi0)) == pval .* pi0

        ## compare with reference values
        @test isapprox( adjust(pval1, t(pi0)), ref, atol = 1e-8 )

        ## BH Adaptive same as BH for π₀ missing or 1
        @test isapprox( adjust(pval1, t(1.0)), ref0, atol = 1e-8 )
        @test isapprox( adjust(pval1, t()), ref0, atol = 1e-8 )
        @test isapprox( adjust(pval1, t(0.0)), zero(ref0), atol = 1e-8 )

        # adaptive with pi0 estimator == oracle with estimated pi0
        @test adjust(pval1, t(Storey())) == adjust(pval1, t(estimate_pi0(pval1, Storey())))
        @test adjust(pval1, t(Storey(0.3))) == adjust(pval1, t(estimate_pi0(pval1, Storey(0.3))))
        @test adjust(pval1, t(LeastSlope())) == adjust(pval1, t(estimate_pi0(pval1, LeastSlope())))

    end

end

end
