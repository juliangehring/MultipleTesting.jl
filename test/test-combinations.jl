### Test combinations ###
module Test_combinations

using MultipleTesting
using StatsBase
using Test


@testset "p-Value combinations" begin

    p1 = [0.01, 0.05, 0.2, 0.8]
    ref1 = Dict(
        Fisher   => 0.01558752, # `metap::sumlog(p)`
        Logit    => 0.020031,   # `metap::logitp(p)`
        Stouffer => 0.02353884, # `metap::sumz(p)`
        Tippett  => 0.03940399, # `gmeta::Cpvaluecombine(p, "tippett")`
        Simes    => 0.04
    )

    p2 = [0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.75]
    ref2 = Dict(
        Fisher   => 1.285885e-06, # `metap::sumlog(p)`
        Logit    => 1.801132e-06, # `metap::logitp(p)`
        Stouffer => 5.092143e-06, # `metap::sumz(p)`
        Tippett  => 0.0007997201, # `gmeta::Cpvaluecombine(p, "tippett")`
        Simes    => 0.0008
    )

    p3 = copy(p1)
    w3 = [0.1, 0.2, 0.3, 0.4] .* 2
    ref3 = Dict(Stouffer => 0.1916892)

    p_adjust_combinations = Dict(
        Sidak             => Tippett,
        BenjaminiHochberg => Simes
    )

    p_single = rand(1)

    # invalid p-value inputs
    p1_invalid = [-1.0, 0.5]
    p2_invalid = [0.5, 1.5]


    @testset "$(method)" for method in keys(ref2)

        @test method <: PValueCombination

        ref = ref1[method]
        @test isapprox(combine(PValues(p1), method()), ref, atol = 1e-8)
        @test isapprox(combine(p1, method()), ref, atol = 1e-8)

        ref = ref2[method]
        @test isapprox(combine(PValues(p2), method()), ref, atol = 1e-8)
        @test isapprox(combine(p2, method()), ref, atol = 1e-8)

        @test_throws DomainError combine(PValues(p1_invalid), method())
        @test_throws DomainError combine(p1_invalid, method())
        @test_throws DomainError combine(PValues(p2_invalid), method())
        @test_throws DomainError combine(p2_invalid, method())

        @test combine(PValues(p_single), method()) == p_single[1]
        @test combine(p_single, method()) == p_single[1]

        for T in (Float32, Float64)
            pv = rand(T, 1)
            @test isa(combine(PValues(pv), method()), T)

            pv = Vector{T}(p1)
            @test isa(combine(PValues(pv), method()), T)
        end

    end


    @testset "Wilkinson combination" begin

        method = Wilkinson(1)

        @test_throws MethodError Wilkinson()
        @test_throws ArgumentError Wilkinson(0)

        @test typeof(method) <: PValueCombination

        # Wilkinson with rank = 1 is Tippett's method
        ref = ref1[Tippett]
        @test isapprox(combine(PValues(p1), method), ref, atol = 1e-8)
        @test isapprox(combine(p1, method), ref, atol = 1e-8)

        ref = ref2[Tippett]
        @test isapprox(combine(PValues(p2), method), ref, atol = 1e-8)
        @test isapprox(combine(p2, method), ref, atol = 1e-8)

        @test_throws ArgumentError combine(PValues(rand(5)), Wilkinson(6))
        @test_throws ArgumentError combine(rand(5), Wilkinson(6))

        @test_throws DomainError combine(PValues(p1_invalid), method)
        @test_throws DomainError combine(p1_invalid, method)
        @test_throws DomainError combine(PValues(p2_invalid), method)
        @test_throws DomainError combine(p2_invalid, method)

        @test combine(PValues(p_single), method) == p_single[1]
        @test combine(p_single, method) == p_single[1]

        for T in (Float32, Float64)
            pv = rand(T, 1)
            @test isa(combine(PValues(pv), method), T)

            pv = Vector{T}(p1)
            @test isa(combine(PValues(pv), method), T)
        end

        # reference values computed with `metap::wilkinsonp(p, r, alpha = 1e-16)$p`
        ref = [0.03940399, 0.01401875, 0.0272, 0.4096]
        p = [combine(PValues(p1), Wilkinson(r)) for r in 1:length(p1)]
        @test isapprox(p, ref, atol = 1e-7)

        ref = [7.997201e-04, 2.788821e-05, 5.393332e-05, 3.717514e-04,
               4.316500e-04, 1.231360e-03, 8.519680e-03, 1.001129e-01]
        p = [combine(p2, Wilkinson(r)) for r in 1:length(p2)]
        @test isapprox(p, ref, atol = 1e-7)

    end


    @testset "Minimum combination with $(p_adjustment)" for (p_adjustment, p_combination) in p_adjust_combinations

        @test_throws MethodError Minimum()

        padj_comb = Minimum(p_adjustment())

        @test typeof(padj_comb) <: PValueCombination

        @test isapprox(combine(PValues(p1), padj_comb), ref1[p_combination], atol = 1e-8)
        @test isapprox(combine(p1, padj_comb), ref1[p_combination], atol = 1e-8)
        @test isapprox(combine(PValues(p2), padj_comb), ref2[p_combination], atol = 1e-8)
        @test isapprox(combine(p2, padj_comb), ref2[p_combination], atol = 1e-8)

        @test_throws DomainError combine(PValues(p1_invalid), padj_comb)
        @test_throws DomainError combine(p1_invalid, padj_comb)
        @test_throws DomainError combine(PValues(p2_invalid), padj_comb)
        @test_throws DomainError combine(p2_invalid, padj_comb)

        @test combine(PValues(p_single), padj_comb) == p_single[1]
        @test combine(p_single, padj_comb) == p_single[1]

        for T in (Float32, Float64)
            pv = rand(T, 1)
            @test isa(combine(PValues(pv), padj_comb), T)

            pv = Vector{T}(p1)
            @test isa(combine(PValues(pv), padj_comb), T)
        end

    end


    @testset "Weighted Stouffer" begin

        method = Stouffer

        ref = ref1[method]
        @test isapprox(combine(p1, fill(1.0, size(p1)), method()), ref, atol = 1e-8)
        @test isapprox(combine(p1, Weights(fill(1.0, size((p1)))), method()), ref, atol = 1e-8)
        @test isapprox(combine(PValues(p1), Weights(fill(1.0, size(p1))), method()), ref, atol = 1e-8)

        ref = ref2[method]
        @test isapprox(combine(p2, fill(1.0, size(p2)), method()), ref, atol = 1e-8)
        @test isapprox(combine(p2, Weights(fill(1.0, size(p2))), method()), ref, atol = 1e-8)
        @test isapprox(combine(PValues(p2), Weights(fill(1.0, size(p2))), method()), ref, atol = 1e-8)

        ref = ref3[method]
        w3norm = w3 ./ sum(w3)
        # unnormalised weights
        @test isapprox(combine(p3, w3, method()), ref, atol = 1e-8)
        @test isapprox(combine(p3, Weights(w3), method()), ref, atol = 1e-8)
        # normalised weights
        @test isapprox(combine(p3, w3norm, method()), ref, atol = 1e-8)
        @test isapprox(combine(p3, Weights(w3norm), method()), ref, atol = 1e-8)

        @test_throws DomainError combine(p1_invalid, fill(1.0, size(p1_invalid)), method())
        @test_throws DomainError combine(p1_invalid, Weights(fill(1.0, size(p1_invalid))), method())
        @test_throws DomainError combine(p2_invalid, fill(1.0, size(p2_invalid)), method())
        @test_throws DomainError combine(p2_invalid, Weights(fill(1.0, size(p2_invalid))), method())

        @test combine(p_single, fill(1.0, size(p_single)), method()) == p_single[1]
        @test combine(p_single, Weights(fill(1.0, size(p_single))), method()) == p_single[1]

    end


    @testset "Edge cases" begin

        p0 = [0.01, 0.05, 0.2, 0.8, 0.0]
        p1 = [0.01, 0.05, 0.2, 0.8, 1.0]

        @test isnan(combine(p0, Fisher()))
        @test isapprox(combine(p1, Fisher()), 0.04198529, atol = 1e-8)

        @test isnan(combine(p0, Logit()))
        @test isnan(combine(p1, Logit()))

        @test isnan(combine(p0, Stouffer()))
        @test isnan(combine(p1, Stouffer()))

        @test isnan(combine(p0, fill(1.0, size(p0)), Stouffer()))
        @test isnan(combine(p1, fill(1.0, size(p1)), Stouffer()))

        @test isapprox(combine(p0, Tippett()), 0.0)
        @test isapprox(combine(p1, Tippett()), 0.04900995, atol = 1e-8)

    end


end

end
