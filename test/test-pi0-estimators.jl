## Test π₀ estimator methods ##
module Test_pi0

using MultipleTesting
using Test
using StatsBase


@testset "π₀ estimators" begin

    ## test case: deterministic
    p0 = collect(0.01:0.01:1)
    p1 = p0 .^ 10
    p = [p0; p1] ## unsorted
    pi0 = length(p0) / length(p)

    lambdas = collect(0.05:0.05:0.95)

    unsort = MultipleTesting.unsort


    @testset "Storey π₀" begin

        @test typeof(Storey()) <: Pi0Estimator
        @test typeof(Storey(0.5)) <: Pi0Estimator
        @test estimate_pi0(p, Storey(0.2)) ≈ 0.6
        @test estimate_pi0(p, Storey(0.0)) ≈ 1.0
        @test estimate_pi0(p, Storey(1.0)) ≈ 1.0

        p_unsort = unsort(p)
        @test !issorted(p_unsort)
        @test estimate_pi0(p_unsort, Storey(0.2)) ≈ 0.6
        @test !issorted(p_unsort)

        @test_throws DomainError Storey(-0.1)
        @test_throws DomainError Storey(1.1)
        #@test_throws DomainError estimate_pi0(p, Storey(1.0)) ## CHCK

    end


    @testset "StoreyBootstrap π₀" begin

        q = 0.1

        @test typeof(StoreyBootstrap()) <: Pi0Estimator

        @test estimate_pi0(p, StoreyBootstrap()) ≈ 0.6
        @test estimate_pi0(p0, StoreyBootstrap()) ≈ 1.0
        @test estimate_pi0(p1, StoreyBootstrap()) ≈ 0.15

        @test estimate_pi0(p, StoreyBootstrap(lambdas, q)) ≈ 0.6
        @test estimate_pi0(p0, StoreyBootstrap(lambdas, q)) ≈ 1.0
        @test estimate_pi0(p1, StoreyBootstrap(lambdas, q)) ≈ 0.15

        ## unsorted lambdas
        lambdas_unsort = unsort(lambdas)
        @test !issorted(lambdas_unsort)
        @test estimate_pi0(p, StoreyBootstrap(lambdas_unsort, q)) ≈ 0.6
        @test !issorted(lambdas_unsort)

        @test_throws DomainError StoreyBootstrap(lambdas, -0.1)
        @test_throws DomainError StoreyBootstrap(lambdas, 1.1)
        @test_throws MethodError StoreyBootstrap(0.5)
        @test_throws MethodError StoreyBootstrap(lambdas)
        @test_throws MethodError StoreyBootstrap(0.5, lambdas)

    end


    @testset "LeastSlope π₀" begin

        @test typeof(LeastSlope()) <: Pi0Estimator

        ## checked against structSSI::pi0.lsl
        @test isapprox( estimate_pi0(p, LeastSlope()), 0.62, atol = 1e-2 )
        @test isapprox( estimate_pi0(p0, LeastSlope()), 1.0, atol = 1e-2 )
        @test isapprox( estimate_pi0(p1, LeastSlope()), 0.16, atol = 1e-2 )

        ## check against internal reference implementation
        @test estimate_pi0(p, LeastSlope()) ≈ MultipleTesting.lsl_pi0_vec(p)
        @test estimate_pi0(p0, LeastSlope()) ≈ MultipleTesting.lsl_pi0_vec(p0)
        @test estimate_pi0(p1, LeastSlope()) ≈ MultipleTesting.lsl_pi0_vec(p1)

        @test_throws MethodError LeastSlope(0.1)

        ## unsorted p-values
        p_unsort = unsort(p)
        @test !issorted(p_unsort)
        @test isapprox( estimate_pi0(p_unsort, LeastSlope()), 0.62, atol = 1e-2 )
        @test !issorted(p_unsort)

        p_unsort = unsort(p)
        @test !issorted(p_unsort)
        @test isapprox( MultipleTesting.lsl_pi0_vec(p_unsort), 0.62, atol = 1e-2 )
        @test !issorted(p_unsort)

    end


    @testset "Oracle π₀" begin

        @test estimate_pi0(p, Oracle(0.5)) == 0.5
        @test estimate_pi0(p0, Oracle(0.6)) == 0.6
        @test estimate_pi0(p1, Oracle()) == 1.0

    end


    @testset "TwoStep π₀" begin

        alpha = 0.05

        @test typeof(TwoStep(alpha)) <: Pi0Estimator

        ## checked against mutoss::TSBKY_pi0_est
        @test estimate_pi0(p, TwoStep()) ≈ 0.665
        @test estimate_pi0(p, TwoStep(alpha)) ≈ 0.665
        @test estimate_pi0(p0, TwoStep(alpha)) ≈ 1.0
        @test estimate_pi0(p1, TwoStep(alpha)) ≈ 0.29
        @test estimate_pi0(p, TwoStep(alpha, BenjaminiHochberg())) ≈ 0.665

        @test estimate_pi0(p, TwoStep(0.1)) ≈ 0.63
        @test estimate_pi0(p, TwoStep(0.0)) ≈ 1.0
        @test estimate_pi0(p, TwoStep(1.0)) ≈ 0.415
        @test estimate_pi0(p, TwoStep(0.1, BenjaminiHochberg())) ≈ 0.63

        ## unsorted p-values
        p_unsort = unsort(p)
        @test !issorted(p_unsort)
        @test estimate_pi0(p_unsort, TwoStep(alpha)) ≈ 0.665
        @test !issorted(p_unsort)

    end


    @testset "RightBoundary π₀" begin

        # R code tested against (note this uses equidistant grid so we have to use
        # λseq as in Storey for comparison
        # ```R
        # library(pi0)
        # p0 <- seq(0.01,1, by=0.01)
        # histf1(p0, rightBoundary=TRUE)
        # [1] 1
        # p1 <- p0^10
        # histf1(p1, rightBoundary=TRUE)
        # [1] 0.1428571
        # p <- c(p0,p1)# histf1(p, rightBoundary=TRUE)
        # [1] 0.5714286
        # ```

        @test typeof(RightBoundary()) <: Pi0Estimator
        @test typeof(RightBoundary(lambdas)) <: Pi0Estimator

        # only eps because pi0 package uses right closed histograms
        @test isapprox( estimate_pi0(p, RightBoundary(lambdas)), 0.5714286, atol = 0.02 )
        @test estimate_pi0(p0, RightBoundary(lambdas)) ≈ 1.0
        @test isapprox( estimate_pi0(p1, RightBoundary(lambdas)), 0.1428571, atol = 1e-7 )

        @test estimate_pi0(p, RightBoundary(lambdas)) ≈ estimate_pi0(p, RightBoundary(unsort(lambdas)))

        # not checked against R implementation but should hold (check for default lambda grid)
        @test estimate_pi0(p0, RightBoundary()) ≈ 1.0

        p_unsort = unsort(p1)
        @test !issorted(p_unsort)
        @test isapprox( estimate_pi0(p_unsort, RightBoundary(lambdas)), 0.1428571, atol = 1e-7 )
        @test !issorted(p_unsort)

    end


    @testset "CensoredBUM π₀" begin

        @test typeof(CensoredBUM()) <: Pi0Estimator
        @test typeof(CensoredBUM(0.2, 0.1)) <: Pi0Estimator

        @test isapprox( estimate_pi0(p, CensoredBUM()), 0.55797, atol = 1e-5 )
        @test isapprox( estimate_pi0(p0, CensoredBUM()), 1.0, atol = 1e-5 )
        @test isapprox( estimate_pi0(p1, CensoredBUM()), 0.11608, atol = 2e-5 )

        # test against internal reference implementation
        # p0 case is not handled well by naive implementation
        @test estimate_pi0(p, CensoredBUM()) ≈ MultipleTesting.cbum_pi0_naive(p)[1]
        @test estimate_pi0(p1, CensoredBUM()) ≈ MultipleTesting.cbum_pi0_naive(p1)[1]

        ## test case that does not converge
        @test isnan(estimate_pi0(p, CensoredBUM(0.5, 0.05, 1e-6, 2)))

        @test_throws DomainError CensoredBUM(-0.5, 0.05)
        @test_throws DomainError CensoredBUM(1.5, 0.05)
        @test_throws DomainError CensoredBUM(0.5, -0.05)
        @test_throws DomainError CensoredBUM(0.5, 1.05)

        @test_throws DomainError CensoredBUM(0.5, 0.05, -1e-6, 100)
        @test_throws DomainError CensoredBUM(0.5, 0.05, 1.1, 100)
        @test_throws DomainError CensoredBUM(0.5, 0.05, 1e-6, -10)

        f = fit(CensoredBUM(), p)
        @test typeof(f) <: CensoredBUMFit
        @test isapprox( f.π0, 0.55797, atol = 1e-5 )

        # denominator becomes 0 if all p-values are 1
        @test isnan( estimate_pi0(ones(20), CensoredBUM()) )

        pi0_est, pars, is_converged = MultipleTesting.cbum_pi0(ones(20))
        @test isnan( pi0_est )
        @test !is_converged

        pi0_est, pars, is_converged = MultipleTesting.cbum_pi0_naive(ones(20))
        @test isnan( pi0_est )
        @test !is_converged

    end


    @testset "BUM π₀" begin

        @test isapprox( estimate_pi0(p, BUM(0.5)), 0.55528, atol = 1e-5 )
        @test isapprox( estimate_pi0(p, BUM()), 0.55528, atol = 1e-5 )
        @test isapprox( estimate_pi0(p0, BUM()), 1.0, atol = 1e-5 )
        @test isapprox( estimate_pi0(p1, BUM()), 0.10874, atol = 1e-5 )

        ## test case that does not converge
        @test isnan(estimate_pi0(p, BUM(0.5, 1e-6, 2)))

        @test typeof(BUM()) <: Pi0Estimator

        @test_throws DomainError BUM(-0.5)
        @test_throws DomainError BUM(1.5)

    end


    @testset "FlatGrenander π₀" begin

        @test typeof(FlatGrenander()) <: Pi0Estimator

        pu = collect(0.1:0.05:0.9)
        @test estimate_pi0(pu, FlatGrenander()) ≈ 1.0
        @test estimate_pi0(pu.^0.5, FlatGrenander()) ≈ 1.0

        # no reference values from literature/software available
        @test estimate_pi0(p0, FlatGrenander()) ≈ 1.0
        @test isapprox( estimate_pi0(p1, FlatGrenander()), 0.10458, atol = 1e-4 )
        @test isapprox( estimate_pi0(p, FlatGrenander()), pi0, atol = 0.1)

        ## longest constant interval: low level
        lci = MultipleTesting.longest_constant_interval

        px = [0:0.1:1;];
        f = [1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.1];
        @test lci(px, f) ≈ 0.5

        f = [1.5, 1.5, 1.5, 1.5, 1.5, 1.2, 1.2, 1.2, 0.2, 0.2, 0.1];
        @test lci(px, f) ≈ 0.2

        f = [0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1];
        @test lci(px, f) ≈ 0.2

        f = [0.5, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1];
        @test lci(px, f) ≈ 0.5

        px = [0.1, 0.3, 0.5, 0.9];
        f = [0.5, 0.5, 0.2, 0.2];
        @test lci(px, f) ≈ 0.2

        px = [0.1, 0.5, 0.7, 0.9];
        f = [0.5, 0.5, 0.2, 0.2];
        @test lci(px, f) ≈ 0.5

    end


    @testset "ConvexDecreasing π₀" begin

        @test typeof(ConvexDecreasing()) <: Pi0Estimator
        @test typeof(ConvexDecreasing(100, 1e-6, 1000)) <: Pi0Estimator

        @test isapprox( estimate_pi0(p, ConvexDecreasing()), 0.5739054, atol = 1e-4)
        @test isapprox( estimate_pi0(p0, ConvexDecreasing()), 1.0, atol = 1e-6 )
        @test isapprox( estimate_pi0(p1, ConvexDecreasing()), 0.1323284, atol = 1e-3)

        @test isapprox( estimate_pi0(ones(20), ConvexDecreasing()), 1.0, atol = 1e-6 )

        ## test case that does not converge
        @test isnan(estimate_pi0(p, ConvexDecreasing(100, 1e-6, 2)))

        @test_throws DomainError ConvexDecreasing(100, 2.0, 1000)
        @test_throws DomainError ConvexDecreasing(100, -1.0, 1000)
        @test_throws DomainError ConvexDecreasing(-100, 0.01, 1000)
        @test_throws DomainError ConvexDecreasing(100, 0.01, -1000)

        f = fit(ConvexDecreasing(), p)
        @test typeof(f) <: ConvexDecreasingFit
        @test isapprox( f.π0, 0.5739054, atol = 1e-4 )

    end

end

end
