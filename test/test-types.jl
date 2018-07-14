## Test types ##
module Test_types

using MultipleTesting
using Test


@testset "Types" begin

    @testset "PValues type" begin

        n = 10
        vals = rand(n)

        pv = PValues(vals)

        # basic vector functionality
        @test values(pv) == vals
        @test sum(pv) == sum(vals)
        @test length(pv) == n
        @test zero(pv) == zeros(eltype(vals), n)
        @test minimum(pv) == minimum(vals)
        @test maximum(pv) == maximum(vals)
        @test extrema(pv) == extrema(vals)
        @test [x for x in pv] == vals
        @test convert(PValues, vals) == PValues(vals)

        # min/max are taken from the type fields,
        # not recomputed from the values
        pvt = PValues(vals)
        @test_throws ErrorException pvt.min = 0.0
        @test_throws ErrorException pvt.max = 1.0

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

        # TODO discuss expected behaviour in edge cases
        @test_throws MethodError PValues(0.5)
        @test_throws ArgumentError PValues([])
        @test_throws TypeError PValues([0, 1])

    end


    @testset "PValues immutability" begin

        n = 10
        vals = rand(n)
        ref = deepcopy(vals)

        pv = PValues(vals)

        @test pv.values == ref
        @test extrema(pv) == extrema(pv.values)

        vals[1] = 0.0
        vals[2] = 1.0

        @test pv.values == ref
        @test pv.values != vals
        @test extrema(pv) == extrema(pv.values)

        @test_throws ErrorException pv[1] = 0.0

    end

end

end
