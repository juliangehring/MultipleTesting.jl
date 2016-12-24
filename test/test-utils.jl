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

        @test_throws DomainError MultipleTesting.validPValues([-1.])
        @test_throws DomainError MultipleTesting.validPValues([2.])

    end


    @testset "reorder" begin

        x = [1, 5, 4, 2, 4, 3]
        no, oo = MultipleTesting.reorder(x)
        @test x[no] == sort(x)
        @test x[no][oo] == x

    end

end

end
