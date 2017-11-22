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
        @test values(pv) == vals
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
        @test_throws TypeError PValues([0, 1])
        @test_throws ArgumentError PValues([]) # Any
        @test_throws ArgumentError PValues(Float64[]) # Float

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


    @testset "ZScores type" begin

        n = 10
        vals = randn(n)

        zs = ZScores(vals)

        # basic vector functionality
        @test values(zs) ≡ vals
        @test sum(zs) == sum(vals)
        @test length(zs) == n
        @test ones(zs) == ones(eltype(vals), n)
        @test minimum(zs) == minimum(vals)
        @test maximum(zs) == maximum(vals)
        @test extrema(zs) == extrema(vals)
        @test [x for x in zs] == vals
        @test convert(ZScores, vals) == ZScores(vals)

        # parametric type that keeps input float type
        for T in (Float16, Float32, Float64)
            vals = randn(T, n)
            zs = ZScores(vals)
            @test eltype(zs) == T
            @test eltype(values(zs)) == T
        end

        # TODO discuss expected behaviour in edge cases
        @test_throws MethodError ZScores(0.5)
        @test_throws TypeError ZScores([0, 1])
        # differences from PValues behaviour
        @test_throws TypeError ZScores([]) # Any
        @test isa( ZScores(Float64[]), ZScores ) # Float

    end


    @testset "PValues - ZScores transformations" begin

        z = [0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0]
        pu = [0.5, 0.15866, 0.84135, 0.02275, 0.97725, 0.00135, 0.99865]

        zs = ZScores(z)
        pv = PValues(pu)

        p1 = PValues(zs, Upper)
        @test typeof(p1) <: PValues
        @test isapprox( p1, pu, atol = 1e-5  )
        z1 = ZScores(p1, Upper)
        @test typeof(z1) <: ZScores
        @test isapprox( z1, z )

        p1 = PValues(zs, Lower)
        @test typeof(p1) <: PValues
        @test isapprox( p1, 1-pu, atol = 1e-5  )
        z1 = ZScores(p1, Lower)
        @test typeof(z1) <: ZScores
        @test isapprox( z1, z )

        p1 = PValues(zs, Both)
        @test typeof(p1) <: PValues
        @test isapprox( p1, 2*min.(pu, 1-pu), atol = 1e-4 )
        # no transformation since sign of z-scores is unknown
        @test_throws MethodError ZScores(p1, Both)

        # defaults and alternatives
        @test PValues(zs, Both) == PValues(zs)
        @test PValues(zs, Both) == PValues(zs, Both())
        @test ZScores(pv, Lower) == ZScores(pv, Lower())

        # only meaningful transformations
        @test_throws MethodError ZScores(zs, Lower)
        @test_throws MethodError PValues(pv, Lower)

        # TODO define desired behaviour
        @test ZScores(zs) == zs
        @test PValues(pv) == pv

    end

end

end
