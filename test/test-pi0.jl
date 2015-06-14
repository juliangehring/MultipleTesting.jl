## Test π0 estimator methods ##
module Test_pi0

using MultipleTesting
using Base.Test
using StatsBase

## test case: deterministic
p0 = [0.01:0.01:1; ]
p1 = p0 .^ 10
p = [p0; p1] ## unsorted
pi0 = length(p0) / length(p)

lambdas = [0.05:0.05:0.95; ]
lambdas_unsort = sample(lambdas, length(lambdas), replace = false)

## storey_pi0 ##
println(" ** ", "storey_pi0")
@test_throws MethodError storey_pi0()
@test_approx_eq storey_pi0(p, 0.2) 0.6
@test_approx_eq storey_pi0(p, 0.) 1.0
#@test_throws DomainError storey_pi0(p, 1.) ## CHCK

## bootstrap_pi0 ##
println(" ** ", "bootstrap_pi0")
@test_throws MethodError bootstrap_pi0()
@test_approx_eq bootstrap_pi0(p, lambdas) 0.6
@test_approx_eq bootstrap_pi0(p0, lambdas) 1.0
@test_approx_eq bootstrap_pi0(p1, lambdas) 0.15
## unsorted lambdas
@test_approx_eq bootstrap_pi0(p, lambdas_unsort) 0.6
## with default 'lambdas'
@test_approx_eq bootstrap_pi0(p) 0.6
@test_approx_eq bootstrap_pi0(p0) 1.0
@test_approx_eq bootstrap_pi0(p1) 0.15

## lsl_pi0 ##
println(" ** ", "lsl_pi0")
@test_approx_eq lsl_pi0(p) MultipleTesting.lsl_pi0_vec(p)
@test_approx_eq lsl_pi0(p0) MultipleTesting.lsl_pi0_vec(p0)
@test_approx_eq lsl_pi0(p1) MultipleTesting.lsl_pi0_vec(p1)

## checked against structSSI::pi0.lsl
@test_approx_eq_eps lsl_pi0(p) 0.62 1e-2
@test_approx_eq_eps lsl_pi0(p0) 1.0 1e-2
@test_approx_eq_eps lsl_pi0(p1) 0.16 1e-2

end
