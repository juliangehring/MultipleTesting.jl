## Test grenander and related methods ##
module Test_grenander

using MultipleTesting
using Base.Test

## Grenander ##
println(" ** ", "grenander")
pv = [0.2, 0.4, 0.6, 0.8, 1.0]
p, f, F = MultipleTesting.grenander(pv)
@test_approx_eq p pv
@test_approx_eq f ones(p)
@test_approx_eq F pv

pv = [0.1, 0.3, 0.5, 0.7, 0.9]
p, f, F = MultipleTesting.grenander(pv)
@test_approx_eq p pv
@test_approx_eq f ones(p)
@test_approx_eq F pv+0.1

pv = [0.1, 0.1, 0.3, 0.3, 0.5, 0.5, 0.7, 0.7, 0.9, 0.9]
p, f, F = MultipleTesting.grenander(pv)
@test_approx_eq p unique(pv)
@test_approx_eq f ones(p)
@test_approx_eq F p+0.1

pv = [0.1, 0.2, 0.5, 0.6, 0.9]
p, f, F = MultipleTesting.grenander(pv)
@test_approx_eq p pv
@test_approx_eq f [2.0, 1.0, 1.0, 2/3, 2/3]
@test_approx_eq F [0.2, 0.4, 0.7, 0.8, 1.0]


## isotonic regression ##
println(" ** ", "isotonic regression")

#pooled adjacent violators example page 10 robertson
isotonicregression = MultipleTesting.isotonicregression
@test_approx_eq isotonicregression([22.5; 23.333; 20.833; 24.25], [3.0;3.0;3.0;2.0]) [22.222; 22.222; 22.222; 24.25]
#if input already ordered, then output should be the same
@test_approx_eq isotonicregression([1.0; 2.0; 3.0]) [1.0; 2.0; 3.0]
@test_approx_eq isotonicregression([1., 41., 51., 1., 2., 5., 24.], [1., 2., 3., 4., 5., 6., 7.]) [1.0, 13.95, 13.95, 13.95, 13.95, 13.95, 24]

#pooled adjacent violators example page 10 robertson
isotonic_regression = MultipleTesting.isotonic_regression
@test_approx_eq isotonic_regression([22.5; 23.333; 20.833; 24.25], [3.0;3.0;3.0;2.0]) [22.222; 22.222; 22.222; 24.25]
#if input already ordered, then output should be the same
@test_approx_eq isotonic_regression([1.0; 2.0; 3.0]) [1.0; 2.0; 3.0]
@test_approx_eq isotonic_regression([1., 41., 51., 1., 2., 5., 24.], [1., 2., 3., 4., 5., 6., 7.]) [1.0, 13.95, 13.95, 13.95, 13.95, 13.95, 24]

# single value or empty vector remains unchanged
r = rand(1)
@test MultipleTesting.isotonic_regression(r) == r
r = Vector{Float64}()
@test MultipleTesting.isotonic_regression(r) == r

r = rand(10)
@test MultipleTesting.isotonic_regression(r, ones(r)) == MultipleTesting.isotonic_regression(r)
@test_throws DimensionMismatch MultipleTesting.isotonic_regression(rand(10), ones(5))


end
