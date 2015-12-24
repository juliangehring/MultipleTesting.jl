## Test pvalue adjustment methods ##
module Test_pval_pi0_adjustment

using MultipleTesting
using Base.Test

pval1 = [0.0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.75, 1.0]
pi0 = 0.4

## compared against:
## qvalue::qvalue(p, pi0, pfdr = TRUE, pi0.method = "bootstrap")

## benjamini_hochberg
m = benjamini_hochberg
ref = [0.0, 0.0005, 0.003333333, 0.025, 0.1, 0.166666667, 0.285714286, 0.5, 0.833333333, 1.0] .* pi0
println(" ** ", m)
@test_throws MethodError m()
## no integers as input
@test_throws MethodError m([0, 1])
## no valid p-values as input
@test_throws DomainError m([-1.0, 0.5], pi0)
@test_throws DomainError m([0.5, 1.5], pi0)
@test_throws DomainError m([0.5, 0.7], -1.0)
@test_throws DomainError m([0.5, 0.7], 1.5)
## single p-value is returned unchanged
pval = rand(1)
@test m(pval, pi0) == pval .* pi0
## compare with reference values
@test_approx_eq_eps m(pval1, pi0) ref 1e-6

## qvalue
m = qValues
ref = [0.0, 0.0002, 0.00133333, 0.01, 0.04, 0.0666667, 0.114286, 0.2, 0.333333, 0.4]
println(" ** ", m)
@test_throws MethodError m()
## no integers as input
@test_throws MethodError m([0, 1])
## no valid p-values as input
@test_throws DomainError m([-1.0, 0.5], pi0)
@test_throws DomainError m([0.5, 1.5], pi0)
@test_throws DomainError m([0.5, 0.7], -1.0)
@test_throws DomainError m([0.5, 0.7], 1.5)
## single p-value is returned unchanged
pval = rand(1)
@test m(pval, pi0) == pval .* pi0
## compare with reference values
@test_approx_eq_eps m(pval1, pi0) ref 1e-6
@test_approx_eq_eps m(pval1, pi0, false) ref 1e-6

## qvalue pfdr
m = qValues
ref = [0.099685, 0.099685, 0.099685, 0.099685, 0.099685, 0.102356, 0.128033, 0.201217, 0.333333, 0.4]
println(" ** ", m)
@test_throws MethodError m()
## no integers as input
@test_throws MethodError m([0, 1])
## no valid p-values as input
@test_throws DomainError m([-1.0, 0.5], pi0)
@test_throws DomainError m([0.5, 1.5], pi0)
@test_throws DomainError m([0.5, 0.7], -1.0)
@test_throws DomainError m([0.5, 0.7], 1.5)
## single p-value is returned unchanged
pval = rand(1)
@test m(pval, pi0) == pval .* pi0
## compare with reference values
@test_approx_eq_eps m(pval1, pi0, true) ref 1e-6

end
