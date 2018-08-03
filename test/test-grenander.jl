## Test grenander and related methods ##
module Test_grenander

using MultipleTesting
using Test


@testset "Grenander estimators" begin

    @testset "Grenander estimator" begin

        pv = [0.2, 0.4, 0.6, 0.8, 1.0]
        p, f, F = MultipleTesting.grenander(pv)
        @test p ≈ pv
        @test f ≈ fill(1.0, size(p))
        @test F ≈ pv

        pv = [0.1, 0.3, 0.5, 0.7, 0.9]
        p, f, F = MultipleTesting.grenander(pv)
        @test p ≈ pv
        @test f ≈ fill(1.0, size(p))
        @test F ≈ pv .+ 0.1

        pv = [0.1, 0.1, 0.3, 0.3, 0.5, 0.5, 0.7, 0.7, 0.9, 0.9]
        p, f, F = MultipleTesting.grenander(pv)
        @test p ≈ unique(pv)
        @test f ≈ fill(1.0, size(p))
        @test F ≈ p .+ 0.1

        pv = [0.1, 0.2, 0.5, 0.6, 0.9]
        p, f, F = MultipleTesting.grenander(pv)
        @test p ≈ pv
        @test f ≈ [2.0, 1.0, 1.0, 2/3, 2/3]
        @test F ≈ [0.2, 0.4, 0.7, 0.8, 1.0]

    end


    @testset "Isotonic regression" begin

        #pooled adjacent violators example page 10 robertson
        isotonic_regression_reference = MultipleTesting.isotonic_regression_reference
        @test isotonic_regression_reference([22.5; 23.333; 20.833; 24.25], [3.0;3.0;3.0;2.0]) ≈ [22.222; 22.222; 22.222; 24.25]

        #if input already ordered, then output should be the same
        @test isotonic_regression_reference([1.0; 2.0; 3.0]) ≈ [1.0; 2.0; 3.0]
        @test isotonic_regression_reference([1., 41., 51., 1., 2., 5., 24.], [1., 2., 3., 4., 5., 6., 7.]) ≈ [1.0, 13.95, 13.95, 13.95, 13.95, 13.95, 24.]

        # single value or empty vector remains unchanged
        r = rand(1)
        @test MultipleTesting.isotonic_regression(r) == r
        r = Vector{Float64}()
        @test MultipleTesting.isotonic_regression(r) == r

        r = rand(10)
        @test MultipleTesting.isotonic_regression(r, fill(1.0, size(r))) == MultipleTesting.isotonic_regression(r)
        @test_throws DimensionMismatch MultipleTesting.isotonic_regression(rand(10), ones(5))

    end

end

end
