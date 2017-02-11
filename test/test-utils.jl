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


    @testset "valid_pvalues" begin

        valid_pvalues = MultipleTesting.valid_pvalues

        # test against invalid vector inputs
        @test_throws DomainError valid_pvalues([-1.])
        @test_throws DomainError valid_pvalues([2.])
        @test_throws DomainError valid_pvalues([0.1, 0.2, 1.2])
        @test_throws DomainError valid_pvalues([0.1, -0.3, 0.2])

        # test against valid vector inputs
        @test valid_pvalues([0.0]) == nothing
        @test valid_pvalues([1.0]) == nothing
        @test valid_pvalues(rand(1)) == nothing
        @test valid_pvalues([0.1, 0.2, 0.9]) == nothing
        @test valid_pvalues(rand(5)) == nothing

    end


    @testset "reorder" begin

        x = [1, 5, 4, 2, 4, 3]
        no, oo = MultipleTesting.reorder(x)
        @test x[no] == sort(x)
        @test x[no][oo] == x

    end


    @testset "sort_if_needed" begin

        x = rand(20)

        sort_if_needed = MultipleTesting.sort_if_needed
        sort_if_needed! = MultipleTesting.sort_if_needed!

        # behaves as standard `sort`
        @test sort_if_needed(x) == sort(x)
        @test sort_if_needed(sort(x)) == sort(x)

        # `sort` keywords work
        @test sort_if_needed(x, rev = true) == sort(x, rev = true)
        @test sort_if_needed(sort(x, rev = true), rev = true) == sort(x, rev = true)

        x = rand(20)
        y = copy(x)
        sort_if_needed!(y)
        @test y == sort(x)
        sort_if_needed!(y)
        @test y == sort(x)

    end


    @testset "unsort" begin

        n = 20
        xs = sort(rand(n))
        xr = reverse(xs) # reverse sorted
        xu = xs[[1:2:n-1; 2:2:n]] # unsorted
        @test !issorted(xu)

        unsort = MultipleTesting.unsort

        @test issorted(xs)
        @test !issorted(unsort(xs))

        # unsorted input gets returned unchanged
        @test unsort(xu) == xu

        # `sort` keywords work
        @test !issorted(xr)
        @test issorted(xr, rev = true)
        @test issorted(unsort(xr), rev = true)
        @test !issorted(unsort(xr, rev = true), rev = true)

    end


end

end
