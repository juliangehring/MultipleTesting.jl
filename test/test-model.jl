## Test models ##
module Test_model

using MultipleTesting
using Test
using Distributions


@testset "Models" begin

    @testset "BUM model" begin

        bum = BetaUniformMixtureModel(0.7)
        bum0 = BetaUniformMixtureModel(1.0) ## uniform, all H₀
        bum1 = BetaUniformMixtureModel(0.0) ## beta, all H₁

        x = collect(0:0.01:1)

        @test cdf.(Ref(bum0), x) ≈ cdf.(Ref(Uniform()), x)
        @test cdf.(Ref(bum1), x) ≈ cdf.(Ref(Beta(0.5, 1.0)), x)

        @test pdf.(Ref(bum0), x) ≈ pdf.(Ref(Uniform()), x)
        @test pdf.(Ref(bum1), x) ≈ pdf.(Ref(Beta(0.5, 1.0)), x)

        @test isin(rand(bum, 100), 0.0, 1.0)
        @test isin(rand(bum0, 100), 0.0, 1.0)
        @test isin(rand(bum1, 100), 0.0, 1.0)

        @test_throws MethodError BetaUniformMixtureModel()
        @test_throws DomainError BetaUniformMixtureModel(-0.2)
        @test_throws DomainError BetaUniformMixtureModel(1.2, 0.2, 2.0)

    end

end

end
