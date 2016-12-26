## Test π0 estimator methods ##
module Test_pi0

using MultipleTesting
using Base.Test
using StatsBase

@testset "π0 estimators" begin

    ## test case: deterministic
    p0 = collect(0.01:0.01:1)
    p1 = p0 .^ 10
    p = [p0; p1] ## unsorted
    pi0 = length(p0) / length(p)

    lambdas = collect(0.05:0.05:0.95)

    function unsort(x)
        y = copy(x)
        while issorted(y)
            sample!(x, y, replace = false)
        end
        return y
    end


    @testset "storey_pi0" begin

        @test_throws MethodError storey_pi0()
        @test storey_pi0(p, 0.2) ≈ 0.6
        @test storey_pi0(p, 0.0) ≈ 1.0
        #@test_throws DomainError storey_pi0(p, 1.) ## CHCK

        p_unsort = unsort(p)
        @test !issorted(p_unsort)
        @test storey_pi0(p_unsort, 0.2) ≈ 0.6
        @test !issorted(p_unsort)

        @test issubtype(typeof(Storey()), Pi0Estimator)
        @test issubtype(typeof(Storey(0.5)), Pi0Estimator)
        @test estimate_pi0(p, Storey(0.2)) ≈ 0.6
        @test estimate_pi0(p, Storey(0.0)) ≈ 1.0
        @test estimate_pi0(p, Storey(1.0)) ≈ 1.0

        @test_throws DomainError Storey(-0.1)
        @test_throws DomainError Storey(1.1)

    end


    @testset "bootstrap_pi0" begin

        @test_throws MethodError bootstrap_pi0()

        @test bootstrap_pi0(p, lambdas) ≈ 0.6
        @test bootstrap_pi0(p0, lambdas) ≈ 1.0
        @test bootstrap_pi0(p1, lambdas) ≈ 0.15

        ## unsorted lambdas
        lambdas_unsort = unsort(lambdas)
        @test !issorted(lambdas_unsort)
        @test bootstrap_pi0(p, lambdas_unsort) ≈ 0.6
        @test !issorted(lambdas_unsort)

        ## with default 'lambdas'
        @test bootstrap_pi0(p) ≈ 0.6
        @test bootstrap_pi0(p0) ≈ 1.0
        @test bootstrap_pi0(p1) ≈ 0.15

        @test issubtype(typeof(StoreyBootstrap()), Pi0Estimator)
        @test estimate_pi0(p, StoreyBootstrap()) ≈ 0.6
        @test estimate_pi0(p0, StoreyBootstrap()) ≈ 1.0
        @test estimate_pi0(p1, StoreyBootstrap()) ≈ 0.15

        @test_throws DomainError StoreyBootstrap(lambdas, -0.1)
        @test_throws DomainError StoreyBootstrap(lambdas, 1.1)
        @test_throws MethodError StoreyBootstrap(0.5)
        @test_throws MethodError StoreyBootstrap(lambdas)
        @test_throws MethodError StoreyBootstrap(0.5, lambdas)

    end


    @testset "lsl_pi0" begin

        @test lsl_pi0(p) ≈ MultipleTesting.lsl_pi0_vec(p)
        @test lsl_pi0(p0) ≈ MultipleTesting.lsl_pi0_vec(p0)
        @test lsl_pi0(p1) ≈ MultipleTesting.lsl_pi0_vec(p1)

        ## checked against structSSI::pi0.lsl
        @test isapprox( lsl_pi0(p), 0.62, atol = 1e-2 )
        @test isapprox( lsl_pi0(p0), 1.0, atol = 1e-2 )
        @test isapprox( lsl_pi0(p1), 0.16, atol = 1e-2 )

        @test issubtype(typeof(LeastSlope()), Pi0Estimator)
        @test isapprox( estimate_pi0(p, LeastSlope()), 0.62, atol = 1e-2 )
        @test isapprox( estimate_pi0(p0, LeastSlope()), 1.0, atol = 1e-2 )
        @test isapprox( estimate_pi0(p1, LeastSlope()), 0.16, atol = 1e-2 )

        @test_throws MethodError LeastSlope(0.1)

        ## unsorted p-values
        p_unsort = unsort(p)
        @test !issorted(p_unsort)
        @test isapprox( lsl_pi0(p_unsort), 0.62, atol = 1e-2 )
        @test !issorted(p_unsort)

        p_unsort = unsort(p)
        @test !issorted(p_unsort)
        @test isapprox( MultipleTesting.lsl_pi0_vec(p_unsort), 0.62, atol = 1e-2 )
        @test !issorted(p_unsort)

    end


    @testset "oracle" begin

        @test estimate_pi0(p, Oracle(0.5)) == 0.5
        @test estimate_pi0(p0, Oracle(0.6)) == 0.6
        @test estimate_pi0(p1, Oracle()) == 1.0

    end


    @testset "twostep_pi0" begin

        alpha = 0.05

        ## checked against mutoss::TSBKY_pi0_est
        @test twostep_pi0(p, alpha) ≈ 0.665
        @test twostep_pi0(p0, alpha) ≈ 1.0
        @test twostep_pi0(p1, alpha) ≈ 0.29

        @test issubtype(typeof(TwoStep(alpha)), Pi0Estimator)
        @test estimate_pi0(p, TwoStep()) ≈ 0.665
        @test estimate_pi0(p, TwoStep(alpha)) ≈ 0.665
        @test estimate_pi0(p0, TwoStep(alpha)) ≈ 1.0
        @test estimate_pi0(p1, TwoStep(alpha)) ≈ 0.29

        @test estimate_pi0(p, TwoStep(0.1)) ≈ 0.63
        @test estimate_pi0(p, TwoStep(0.0)) ≈ 1.0
        @test estimate_pi0(p, TwoStep(1.0)) ≈ 0.415

        ## unsorted p-values
        p_unsort = unsort(p)
        @test !issorted(p_unsort)
        @test twostep_pi0(p_unsort, alpha) ≈ 0.665
        @test !issorted(p_unsort)

    end


    @testset "rightboundary_pi0" begin

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

        # only eps because pi0 package uses right closed histograms
        @test isapprox( rightboundary_pi0(p, lambdas), 0.5714286, atol = 0.02 )
        @test rightboundary_pi0(p0, lambdas) ≈ 1.0
        @test isapprox( rightboundary_pi0(p1, lambdas), 0.1428571, atol = 1e-7 )

        @test rightboundary_pi0(p, lambdas) ≈ rightboundary_pi0(p, unsort(lambdas))

        @test issubtype(typeof(RightBoundary(lambdas)), Pi0Estimator)
        @test isapprox( estimate_pi0(p, RightBoundary(lambdas)), 0.5714286, atol = 0.02 )
        @test estimate_pi0(p0, RightBoundary(lambdas)) ≈ 1.0
        @test isapprox( estimate_pi0(p1, RightBoundary(lambdas)), 0.1428571, atol = 1e-7 )

        # not checked against R implementation but should hold (check for default lambda grid)
        @test estimate_pi0(p0, RightBoundary()) ≈ 1.0

        @test issubtype(typeof(RightBoundary()), Pi0Estimator)

        p_unsort = unsort(p1)
        @test !issorted(p_unsort)
        @test isapprox( rightboundary_pi0(p_unsort, lambdas), 0.1428571, atol = 1e-7 )
        @test !issorted(p_unsort)

    end


    @testset "censoredBUM_pi0" begin

        @test estimate_pi0(p, CensoredBUM()) ≈ MultipleTesting.cbum_pi0_naive(p)[1]

        @test isapprox( estimate_pi0(p, CensoredBUM()), 0.55797, atol = 1e-5 )
        @test isapprox( estimate_pi0(p0, CensoredBUM()), 1.0, atol = 1e-5 )
        @test isapprox( estimate_pi0(p1, CensoredBUM()), 0.11608, atol = 2e-5 )
        @test isapprox( estimate_pi0(ones(50), CensoredBUM()), 1.0, atol = 1e-5 )

        ## test case that does not converge
        stderr_dump = redirect_stderr()
        @test isnan(estimate_pi0(p, CensoredBUM(0.5, 0.05, 1e-6, 2)))

        @test issubtype(typeof(CensoredBUM()), Pi0Estimator)
        @test issubtype(typeof(CensoredBUM(0.2, 0.1)), Pi0Estimator)

        @test_throws DomainError CensoredBUM(-0.5, 0.05)
        @test_throws DomainError CensoredBUM(1.5, 0.05)
        @test_throws DomainError CensoredBUM(0.5, -0.05)
        @test_throws DomainError CensoredBUM(0.5, 1.05)

        @test_throws DomainError CensoredBUM(0.5, 0.05, -1e-6, 100)
        @test_throws DomainError CensoredBUM(0.5, 0.05, 1.1, 100)
        @test_throws DomainError CensoredBUM(0.5, 0.05, 1e-6, -10)

        f = fit(CensoredBUM(), p)
        @test issubtype(typeof(f), CensoredBUMFit)
        @test isapprox( f.π0, 0.55797, atol = 1e-5 )

        pi0_est, pars, is_converged = MultipleTesting.cbum_pi0(ones(50))
        @test pi0_est ≈ 1.0
        @test is_converged

        pi0_est, pars, is_converged = MultipleTesting.cbum_pi0_naive(ones(50))
        @test isnan(pi0_est)
        @test !is_converged

    end


    @testset "BUM_pi0" begin

        @test isapprox( estimate_pi0(p, BUM(0.5)), 0.55528, atol = 1e-5 )
        @test isapprox( estimate_pi0(p, BUM()), 0.55528, atol = 1e-5 )
        @test isapprox( estimate_pi0(p0, BUM()), 1.0, atol = 1e-5 )
        @test isapprox( estimate_pi0(p1, BUM()), 0.10874, atol = 1e-5 )

        ## test case that does not converge
        @test isnan(estimate_pi0(p, BUM(0.5, 1e-6, 2)))

        @test issubtype(typeof(BUM()), Pi0Estimator)

        @test_throws DomainError BUM(-0.5)
        @test_throws DomainError BUM(1.5)

    end


    @testset "Flat Grenander" begin

        FlatGrenander = MultipleTesting.FlatGrenander

        @test issubtype(typeof(FlatGrenander()), Pi0Estimator)

        pu = collect(0.1:0.05:0.9)
        @test estimate_pi0(pu, FlatGrenander()) ≈ 1.0
        @test estimate_pi0(pu.^0.5, FlatGrenander()) ≈ 1.0

        @test estimate_pi0(p0, FlatGrenander()) ≈ 1.0
        @test estimate_pi0(p1, FlatGrenander()) < 0.15
        @test isapprox( estimate_pi0(p, FlatGrenander()), pi0, atol = 0.1)

        ## longest constant interval: low level
        lci = MultipleTesting.longest_constant_interval

        p = [0:0.1:1;];
        f = [1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.1];
        @test lci(p, f) ≈ 0.5

        f = [1.5, 1.5, 1.5, 1.5, 1.5, 1.2, 1.2, 1.2, 0.2, 0.2, 0.1];
        @test lci(p, f) ≈ 0.2

        f = [0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1];
        @test lci(p, f) ≈ 0.2

        f = [0.5, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1];
        @test lci(p, f) ≈ 0.5

        p = [0.1, 0.3, 0.5, 0.9];
        f = [0.5, 0.5, 0.2, 0.2];
        @test lci(p, f) ≈ 0.2

        p = [0.1, 0.5, 0.7, 0.9];
        f = [0.5, 0.5, 0.2, 0.2];
        @test lci(p, f) ≈ 0.5

    end

end

end
