## Test π0 estimator methods ##
module Test_pi0

using MultipleTesting
using Base.Test
using StatsBase

## test case: deterministic
p0 = collect(0.01:0.01:1)
p1 = p0 .^ 10
p = [p0; p1] ## unsorted
pi0 = length(p0) / length(p)

function unsort(x)
    y = copy(x)
    while issorted(y)
        sample!(x, y, replace = false)
    end
    return y
end

lambdas = collect(0.05:0.05:0.95)

## storey_pi0 ##
println(" ** ", "storey_pi0")
@test_throws MethodError storey_pi0()
@test_approx_eq storey_pi0(p, 0.2) 0.6
@test_approx_eq storey_pi0(p, 0.0) 1.0
#@test_throws DomainError storey_pi0(p, 1.) ## CHCK

p_unsort = unsort(p)
@test !issorted(p_unsort)
@test_approx_eq storey_pi0(p_unsort, 0.2) 0.6
@test !issorted(p_unsort)

@test issubtype(typeof(Storey()), Pi0Estimator)
@test issubtype(typeof(Storey(0.5)), Pi0Estimator)
@test_approx_eq estimate_pi0(p, Storey(0.2)) 0.6
@test_approx_eq estimate_pi0(p, Storey(0.0)) 1.0
@test_approx_eq estimate_pi0(p, Storey(1.0)) 1.0

@test_throws DomainError Storey(-0.1)
@test_throws DomainError Storey(1.1)



## bootstrap_pi0 ##
println(" ** ", "bootstrap_pi0")
@test_throws MethodError bootstrap_pi0()

@test_approx_eq bootstrap_pi0(p, lambdas) 0.6
@test_approx_eq bootstrap_pi0(p0, lambdas) 1.0
@test_approx_eq bootstrap_pi0(p1, lambdas) 0.15

## unsorted lambdas
lambdas_unsort = unsort(lambdas)
@test !issorted(lambdas_unsort)
@test_approx_eq bootstrap_pi0(p, lambdas_unsort) 0.6
@test !issorted(lambdas_unsort)

## with default 'lambdas'
@test_approx_eq bootstrap_pi0(p) 0.6
@test_approx_eq bootstrap_pi0(p0) 1.0
@test_approx_eq bootstrap_pi0(p1) 0.15

@test issubtype(typeof(StoreyBootstrap()), Pi0Estimator)
@test_approx_eq estimate_pi0(p, StoreyBootstrap()) 0.6
@test_approx_eq estimate_pi0(p0, StoreyBootstrap()) 1.0
@test_approx_eq estimate_pi0(p1, StoreyBootstrap()) 0.15

@test_throws DomainError StoreyBootstrap(lambdas, -0.1)
@test_throws DomainError StoreyBootstrap(lambdas, 1.1)
@test_throws MethodError StoreyBootstrap(0.5)
@test_throws MethodError StoreyBootstrap(lambdas)
@test_throws MethodError StoreyBootstrap(0.5, lambdas)

## lsl_pi0 ##
println(" ** ", "lsl_pi0")
@test_approx_eq lsl_pi0(p) MultipleTesting.lsl_pi0_vec(p)
@test_approx_eq lsl_pi0(p0) MultipleTesting.lsl_pi0_vec(p0)
@test_approx_eq lsl_pi0(p1) MultipleTesting.lsl_pi0_vec(p1)

## checked against structSSI::pi0.lsl
@test_approx_eq_eps lsl_pi0(p) 0.62 1e-2
@test_approx_eq_eps lsl_pi0(p0) 1.0 1e-2
@test_approx_eq_eps lsl_pi0(p1) 0.16 1e-2

@test issubtype(typeof(LeastSlope()), Pi0Estimator)
@test_approx_eq_eps estimate_pi0(p, LeastSlope()) 0.62 1e-2
@test_approx_eq_eps estimate_pi0(p0, LeastSlope()) 1.0 1e-2
@test_approx_eq_eps estimate_pi0(p1, LeastSlope()) 0.16 1e-2

@test_throws MethodError LeastSlope(0.1)

## unsorted p-values
p_unsort = unsort(p)
@test !issorted(p_unsort)
@test_approx_eq_eps lsl_pi0(p_unsort) 0.62 1e-2
@test !issorted(p_unsort)

p_unsort = unsort(p)
@test !issorted(p_unsort)
@test_approx_eq_eps MultipleTesting.lsl_pi0_vec(p_unsort) 0.62 1e-2
@test !issorted(p_unsort)

## Oracle ##
println(" ** ", "oracle")
@test estimate_pi0(p, Oracle(0.5)) == 0.5
@test estimate_pi0(p0, Oracle(0.6)) == 0.6
@test estimate_pi0(p1, Oracle()) == 1.0

## twostep_pi0 ##
println(" ** ", "twostep_pi0")

alpha = 0.05
alpha2 = 0.1

## checked against mutoss::TSBKY_pi0_est
@test_approx_eq twostep_pi0(p, alpha) 0.665
@test_approx_eq twostep_pi0(p0, alpha) 1.0
@test_approx_eq twostep_pi0(p1, alpha) 0.29

@test issubtype(typeof(TwoStep(alpha)), Pi0Estimator)
@test_approx_eq estimate_pi0(p, TwoStep(alpha)) 0.665
@test_approx_eq estimate_pi0(p0, TwoStep(alpha)) 1.0
@test_approx_eq estimate_pi0(p1, TwoStep(alpha)) 0.29

@test_approx_eq estimate_pi0(p, TwoStep(0.1)) 0.63
@test_approx_eq estimate_pi0(p, TwoStep(0.0)) 1.0
@test_approx_eq estimate_pi0(p, TwoStep(1.0)) 0.415

@test_throws MethodError TwoStep()

## unsorted p-values
p_unsort = unsort(p)
@test !issorted(p_unsort)
@test_approx_eq twostep_pi0(p_unsort, alpha) 0.665
@test !issorted(p_unsort)

end
