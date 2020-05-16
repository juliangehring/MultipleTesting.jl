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
        @test f ≈ [2.0, 1.0, 1.0, 2 / 3, 2 / 3]
        @test F ≈ [0.2, 0.4, 0.7, 0.8, 1.0]

    end


    @testset "Isotonic regression" begin

        function isotonic_regression_reference(y::AbstractVector{T}, w::AbstractVector{T}) where T <: AbstractFloat
            y = copy(y)
            w = copy(w)
            cnts = fill(1, length(y))
            i = 2
            while i <= length(y)
                if y[i] < y[i - 1]
                    y[i - 1] = (w[i] * y[i] + w[i - 1] * y[i - 1]) / (w[i] + w[i - 1])
                    w[i - 1] = w[i] + w[i - 1]
                    cnts[i - 1] += cnts[i]
                    deleteat!(y, i)
                    deleteat!(w, i)
                    deleteat!(cnts, i)
                    i = max(i - 2, 1)
                end
                i += 1
            end
            yisotonic = vcat([y[idx] .* fill(1.0, cnt) for (idx, cnt) in enumerate(cnts)]...)
            return yisotonic
        end

        function isotonic_regression_reference(y::AbstractVector{T}) where T <: AbstractFloat
            isotonic_regression_reference(y, fill(1.0, size(y)))
        end


        # pooled adjacent violators example: Robertson, page 10
        @test isotonic_regression_reference([22.5; 23.333; 20.833; 24.25], [3.0;3.0;3.0;2.0]) ≈ [22.222; 22.222; 22.222; 24.25]

        # if input already ordered, then output should be the same
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
