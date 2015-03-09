## Test π0 estimator methods ##

using MultipleTesting
using Base.Test

## test case: deterministic
p0 = [0.01:0.01:1; ]
p1 = p0 .^ 10
p = [p0; p1] ## unsorted
pi0 = length(p0) / length(p)

lambdas = [0.05:0.05:0.95; ]

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
