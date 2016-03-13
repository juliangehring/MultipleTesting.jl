## Test util methods ##
module Test_utils

using MultipleTesting
using Base.Test

## isin ##
println(" ** ", "isin")
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

## validPValues ##
println(" ** ", "validPValues")
@test_throws DomainError MultipleTesting.validPValues([-1.])
@test_throws DomainError MultipleTesting.validPValues([2.])

## reorder ##
println(" ** ", "reorder")
x = [1, 5, 4, 2, 4, 3]
no, oo = MultipleTesting.reorder(x)
@test x[no] == sort(x)
@test x[no][oo] == x

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
println(" ** ", "isotonic_regression")
r = rand(10)
@test MultipleTesting.isotonic_regression(r, ones(r)) == MultipleTesting.isotonic_regression(r)
@test_throws DimensionMismatch MultipleTesting.isotonic_regression(rand(10), ones(5))

end
