### Test combinations ###
module Test_combinations

using MultipleTesting
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

    @testset "$(method)" for method in keys(ref2)

        @test issubtype(method, PValueCombinationMethod)
        @test issubtype(typeof(method()), PValueCombinationMethod)

        ref = ref1[method]
        @test isapprox( combine(p1, method()), ref, atol = 1e-8)

        ref = ref2[method]
        @test isapprox( combine(p2, method()), ref, atol = 1e-8)

    end


    @testset "Stouffer combination weighted" begin

        method = StoufferCombination

        ref = ref1[method]
        @test isapprox( combine(p1, ones(p1), method()), ref, atol = 1e-8)

        ref = ref2[method]
        @test isapprox( combine(p2, ones(p2), method()), ref, atol = 1e-8)

        ref = ref3[method]
        w3norm = w3 ./ sum(w3)
        # unnormalised weights
        @test isapprox( combine(p3, w3, method()), ref, atol = 1e-8 )
        # normalised weights
        @test isapprox( combine(p3, w3norm, method()), ref, atol = 1e-8 )

    end

end

end
