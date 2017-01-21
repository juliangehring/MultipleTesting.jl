## Test types ##
module Test_types

using MultipleTesting
using Base.Test

@testset "Types" begin

    @testset "PValues type" begin

        n = 10
        vals = rand(n)

        pv = PValues(vals)

        # basic vector functionality
        @test values(pv) ≡ vals
        @test sum(pv) == sum(vals)
        @test length(pv) == n
        @test ones(pv) == ones(eltype(vals), n)
        @test minimum(pv) == minimum(vals)
        @test maximum(pv) == maximum(vals)
        @test extrema(pv) == extrema(vals)
        @test [x for x in pv] == vals
        @test convert(PValues, vals) == PValues(vals)

        # min/max are taken from the type fields,
        # not recomputed from the values
        pvt = PValues(vals)
        pvt.min = 0.0
        pvt.max = 1.0
        @test minimum(pvt) == 0.0
        @test maximum(pvt) == 1.0

        # parametric type that keeps input float type
        for T in (Float16, Float32, Float64)
            vals = rand(T, n)
            pv = PValues(vals)
            @test eltype(pv) == T
            @test eltype(values(pv)) == T
        end

        # min/max cannot be passed to constructor
        @test_throws MethodError PValues(vals, 0.0, 1.0)

        # invalid p-values are caught
        @test_throws DomainError PValues([0.1, 0.2, 2.0])
        @test_throws DomainError PValues([0.1, 0.2, -2.0])

    end

end

end
