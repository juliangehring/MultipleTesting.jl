## Test models ##
module Test_model

using MultipleTesting
using Base.Test
using Distributions

@testset "Models" begin

    @testset "BUM model" begin

        bum = BetaUniformMixtureModel(0.7)
        bum0 = BetaUniformMixtureModel(1.0) ## uniform, all H_0
        bum1 = BetaUniformMixtureModel(0.0) ## beta, all H_1

        x = collect(0:0.01:1)

        @test_approx_eq cdf(bum0, x) cdf(Uniform(), x)
        @test_approx_eq cdf(bum1, x) cdf(Beta(0.5, 1.0), x)

        @test_approx_eq pdf(bum0, x) pdf(Uniform(), x)
        @test_approx_eq pdf(bum1, x) pdf(Beta(0.5, 1.0), x)

        @test isin(rand(bum, 100), 0.0, 1.0)
        @test isin(rand(bum0, 100), 0.0, 1.0)
        @test isin(rand(bum1, 100), 0.0, 1.0)

        @test_throws MethodError BetaUniformMixtureModel()
        @test_throws DomainError BetaUniformMixtureModel(-0.2)
        @test_throws DomainError BetaUniformMixtureModel(1.2, 0.2, 2.0)

    end

end

end
