## Test util methods ##
module Test_utils

using MultipleTesting
using Base.Test

@testset "Utility functions" begin

    @testset  "isin" begin

        @test isin(0.1)
        @test isin(0.1, 0.0, 1.0)
        @test !isin(2.1)
        @test !isin(2.1, 0.0, 1.0)
        @test isin(rand(10))
        @test isin(rand(10), 0.0, 1.0)
        @test !isin(-rand(10))
        @test !isin(-rand(10), 0.0, 1.0)
        @test !isin([rand(10); -0.1])
        @test !isin([rand(10); -0.1], 0.0, 1.0)

    end


    @testset "validPValues" begin

        validPValues = MultipleTesting.validPValues

        # test against invalid vector inputs
        @test_throws DomainError validPValues([-1.])
        @test_throws DomainError validPValues([2.])
        @test_throws DomainError validPValues([0.1, 0.2, 1.2])
        @test_throws DomainError validPValues([0.1, -0.3, 0.2])

        # test against valid vector inputs
        @test validPValues([0.0]) == nothing
        @test validPValues([1.0]) == nothing
        @test validPValues(rand(1)) == nothing
        @test validPValues([0.1, 0.2, 0.9]) == nothing
        @test validPValues(rand(5)) == nothing

    end


    @testset "reorder" begin

        x = [1, 5, 4, 2, 4, 3]
        no, oo = MultipleTesting.reorder(x)
        @test x[no] == sort(x)
        @test x[no][oo] == x

    end

end

end
