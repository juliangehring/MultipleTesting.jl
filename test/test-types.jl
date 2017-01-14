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
        @test values(pv) === vals
        @test sum(pv) == sum(vals)
        @test length(pv) == n
        @test ones(pv) == ones(eltype(vals), n)
        @test minimum(pv) == minimum(vals)
        @test maximum(pv) == maximum(vals)
        @test extrema(pv) == extrema(vals)
        @test [x for x in pv] == vals

        pvt = PValues(vals, 0.0, 1.0)
        @test minimum(pvt) == 0.0
        @test maximum(pvt) == 1.0

        for T in (Float16, Float32, Float64)
            vals = rand(T, n)
            pv = PValues(vals)
            @test eltype(pv) == T
            @test eltype(values(pv)) == T
        end

        @test_throws DomainError PValues([0.1, 0.2, 2.0])
        @test_throws DomainError PValues([0.1, 0.2, -2.0])

    end

end

end
