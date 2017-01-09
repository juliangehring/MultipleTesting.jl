### Test combinations ###
module Test_combinations

using MultipleTesting
using StatsBase
using Base.Test

@testset "p-Value combinations" begin

    p1 = [0.01, 0.05, 0.2, 0.8]
    ref1 = Dict(
        FisherCombination   => 0.01558752, # metap::sumlog(p)
        LogitCombination    => 0.020031,   # metap::logitp(p)
        StoufferCombination => 0.02353884, # metap::sumz(p)
        TippettCombination  => 0.03940399  # gmeta::Cpvaluecombine(p, "tippett")
    )

    p2 = [0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.75]
    ref2 = Dict(
        FisherCombination   => 1.285885e-06, # metap::sumlog(p)
        LogitCombination    => 1.801132e-06, # metap::logitp(p)
        StoufferCombination => 5.092143e-06, # metap::sumz(p)
        TippettCombination  => 0.0007997201  # gmeta::Cpvaluecombine(p, "tippett")
    )

    p3 = copy(p1)
    w3 = [0.1, 0.2, 0.3, 0.4] .* 2
    ref3 = Dict(
        StoufferCombination => 0.1916892 # metap::sumz(p, w)
    )

    p_single = rand(1)

    # invalid p-value inputs
    p1_invalid = [-1.0, 0.5]
    p2_invalid = [0.5, 1.5]


    @testset "$(method)" for method in keys(ref2)

        @test issubtype(method, PValueCombinationMethod)
        @test issubtype(typeof(method()), PValueCombinationMethod)

        ref = ref1[method]
        @test isapprox( combine(p1, method()), ref, atol = 1e-8)

        ref = ref2[method]
        @test isapprox( combine(p2, method()), ref, atol = 1e-8)

        @test_throws DomainError combine(p1_invalid, method())
        @test_throws DomainError combine(p2_invalid, method())

        @test combine(p_single, method()) == p_single[1]

    end


    @testset "Weighted StoufferCombination" begin

        method = StoufferCombination

        # `WeightVec` and `weights` are the same
        @test WeightVec(w3) ≡ weights(w3)

        ref = ref1[method]
        @test isapprox( combine(p1, ones(p1), method()), ref, atol = 1e-8)
        @test isapprox( combine(p1, weights(ones(p1)), method()), ref, atol = 1e-8)

        ref = ref2[method]
        @test isapprox( combine(p2, ones(p2), method()), ref, atol = 1e-8)
        @test isapprox( combine(p2, weights(ones(p2)), method()), ref, atol = 1e-8)

        ref = ref3[method]
        w3norm = w3 ./ sum(w3)
        # unnormalised weights
        @test isapprox( combine(p3, w3, method()), ref, atol = 1e-8 )
        @test isapprox( combine(p3, weights(w3), method()), ref, atol = 1e-8 )
        # normalised weights
        @test isapprox( combine(p3, w3norm, method()), ref, atol = 1e-8 )
        @test isapprox( combine(p3, weights(w3norm), method()), ref, atol = 1e-8 )

        @test_throws DomainError combine(p1_invalid, ones(p1_invalid), method())
        @test_throws DomainError combine(p1_invalid, weights(ones(p1_invalid)), method())
        @test_throws DomainError combine(p2_invalid, ones(p2_invalid), method())
        @test_throws DomainError combine(p2_invalid, weights(ones(p2_invalid)), method())

        @test combine(p_single, ones(p_single), method()) == p_single[1]
        @test combine(p_single, weights(ones(p_single)), method()) == p_single[1]

    end


    @testset "Edge cases" begin

        p0 = [0.01, 0.05, 0.2, 0.8, 0.0]
        p1 = [0.01, 0.05, 0.2, 0.8, 1.0]

        @test isnan( combine(p0, FisherCombination()) )
        @test isapprox( combine(p1, FisherCombination()), 0.04198529, atol = 1e-8 )

        @test isnan( combine(p0, LogitCombination()) )
        @test isnan( combine(p1, LogitCombination()) )

        @test isnan( combine(p0, StoufferCombination()) )
        @test isnan( combine(p1, StoufferCombination()) )

        @test isnan( combine(p0, ones(p0), StoufferCombination()) )
        @test isnan( combine(p1, ones(p1), StoufferCombination()) )

        @test isapprox( combine(p0, TippettCombination()), 0.0 )
        @test isapprox( combine(p1, TippettCombination()), 0.04900995, atol = 1e-8 )

    end


end

end
