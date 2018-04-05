## Test pvalue adjustment methods ##
module Test_pval_pi0_adjustment

using MultipleTesting
using Base.Test


@testset "p-Value π0 adjustment" begin

    pval1 = [0.0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.75, 1.0]
    pi0 = 0.4

    @test_throws MethodError adjust()
    @test_throws MethodError adjust(pval1)

    @testset "Benjamini-Hochberg" begin

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

        ## BH Adaptive same as BH for π0 missing or 1
        @test isapprox( adjust(pval1, t(1.0)), ref0, atol = 1e-8 )
        @test isapprox( adjust(pval1, t()), ref0, atol = 1e-8 )
        @test isapprox( adjust(pval1, t(0.0)), zeros(ref0), atol = 1e-8 )

        # adaptive with pi0 estimator == oracle with estimated pi0
        @test adjust(pval1, t(Storey())) == adjust(pval1, t(estimate_pi0(pval1, Storey())))
        @test adjust(pval1, t(Storey(0.3))) == adjust(pval1, t(estimate_pi0(pval1, Storey(0.3))))
        @test adjust(pval1, t(LeastSlope())) == adjust(pval1, t(estimate_pi0(pval1, LeastSlope())))

    end


    @testset "q-value" begin

        # default constructor
        @test QValues() == QValues(Oracle(1.0), false)

        # reference values computed with the `qvalue` package, v2.10.0
        # ```qvalue::qvalue(pval1, pi0=pi0, pfdr = FALSE, lfdr.out=FALSE)$qvalues```
        ref = [0.0, 0.0002, 0.001333333, 0.01, 0.04, 0.066666667, 0.114285714, 0.2, 0.333333333, 0.4]

        qvt = QValues(Oracle(pi0), false)

        @test typeof(qvt.pi0estimator) <: Pi0Estimator
        @test qvt.pfdr == false

        ## single p-value is returned unchanged
        for i in 1:10
            pval = PValues(rand(1))
            @test estimate(pval, qvt) == pval .* pi0
        end

        ## compare with reference values
        @test isapprox( estimate(PValues(pval1), qvt), ref, atol = 1e-8 )

        ## test with all π0 estimators, with different p-values
        for pi0est in subtypes(Pi0Estimator)
            for pv in (PValues(pval1), PValues(rand(100)))
                if isin( estimate_pi0(pv, pi0est()) )
                    qvals = estimate(pv, QValues(pi0est(), false))
                    @test isin(qvals, 0.0, 1.0)
                end
            end
        end

    end


    @testset "q-value with pFDR" begin

        # reference values computed with the `qvalue` package, v2.10.0
        # ```qvalue::qvalue(pval1, pi0=pi0, pfdr = TRUE, lfdr.out=FALSE)$qvalues```
        ref = [NaN, 0.09968523, 0.09968523, 0.09968523, 0.09968523, 0.10235600, 0.12803317, 0.20121668, 0.33333365, 0.4]

        qvt = QValues(Oracle(pi0), true)

        @test typeof(qvt.pi0estimator) <: Pi0Estimator
        @test qvt.pfdr == true

        ## single p-value yields a q-value of π0
        for i in 1:10
            pval = PValues(rand(1))
            @test estimate(pval, qvt) ==  [pi0]
        end

        ## compare with reference values
        @test isapprox( estimate(PValues(pval1), qvt), ref, atol = 1e-8, nans = true )

        ## test with all π0 estimators, with different p-values
        for pi0est in subtypes(Pi0Estimator)
            for pv in (PValues(rand(10)), PValues(rand(100)))
                if isin( estimate_pi0(pv, pi0est()) )
                    qvals = estimate(pv, QValues(pi0est(), true))
                    @test isin(qvals, 0.0, 1.0)
                end
            end
        end

    end

end

end
