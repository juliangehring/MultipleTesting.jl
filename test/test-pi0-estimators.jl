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

## checked against mutoss::TSBKY_pi0_est
@test_approx_eq twostep_pi0(p, alpha) 0.665
@test_approx_eq twostep_pi0(p0, alpha) 1.0
@test_approx_eq twostep_pi0(p1, alpha) 0.29

@test issubtype(typeof(TwoStep(alpha)), Pi0Estimator)
@test_approx_eq estimate_pi0(p, TwoStep()) 0.665
@test_approx_eq estimate_pi0(p, TwoStep(alpha)) 0.665
@test_approx_eq estimate_pi0(p0, TwoStep(alpha)) 1.0
@test_approx_eq estimate_pi0(p1, TwoStep(alpha)) 0.29

@test_approx_eq estimate_pi0(p, TwoStep(0.1)) 0.63
@test_approx_eq estimate_pi0(p, TwoStep(0.0)) 1.0
@test_approx_eq estimate_pi0(p, TwoStep(1.0)) 0.415

## unsorted p-values
p_unsort = unsort(p)
@test !issorted(p_unsort)
@test_approx_eq twostep_pi0(p_unsort, alpha) 0.665
@test !issorted(p_unsort)


## rightboundary_pi0 ##
println(" ** ", "rightboundary_pi0")

# R code tested against (note this uses equidistant grid so we have to use
# λseq as in Storey for comparison

# library(pi0)
# p0 <- seq(0.01,1, by=0.01)
# histf1(p0, rightBoundary=TRUE)
# [1] 1
# p1 <- p0^10
# histf1(p1, rightBoundary=TRUE)
# [1] 0.1428571
# p <- c(p0,p1)# histf1(p, rightBoundary=TRUE)
# [1] 0.5714286

# only eps because pi0 package uses right closed histograms
@test_approx_eq_eps rightboundary_pi0(p, lambdas) 0.5714286 0.02
@test_approx_eq rightboundary_pi0(p0, lambdas) 1.0
@test_approx_eq_eps rightboundary_pi0(p1, lambdas) 0.1428571 10.0^(-7)

@test_approx_eq rightboundary_pi0(p, lambdas) rightboundary_pi0(p, unsort(lambdas))

@test issubtype(typeof(RightBoundary(lambdas)), Pi0Estimator)
@test_approx_eq_eps estimate_pi0(p, RightBoundary(lambdas)) 0.5714286 0.02
@test_approx_eq estimate_pi0(p0, RightBoundary(lambdas)) 1.0
@test_approx_eq_eps estimate_pi0(p1, RightBoundary(lambdas)) 0.1428571 10.0^(-7)

# not checked against R implementation but should hold (check for default lambda grid)
@test_approx_eq estimate_pi0(p0, RightBoundary()) 1.0

@test issubtype(typeof(RightBoundary()), Pi0Estimator)

p_unsort = unsort(p1)
@test !issorted(p_unsort)
@test_approx_eq_eps rightboundary_pi0(p_unsort, lambdas) 0.1428571 10.0^(-7)
@test !issorted(p_unsort)


## censoredBUM_pi0 ##
println(" ** ", "censoredBUM_pi0")

@test_approx_eq_eps estimate_pi0(p, CensoredBUM()) 0.55797 1e-5
@test_approx_eq estimate_pi0(p, CensoredBUM()) MultipleTesting.cbum_pi0_naive(p)[1]
@test_approx_eq_eps estimate_pi0(p0, CensoredBUM()) 1.0 1e-5
@test_approx_eq_eps estimate_pi0(p1, CensoredBUM()) 0.11608 2e-5
@test_approx_eq_eps estimate_pi0(ones(50), CensoredBUM()) 1.0 1e-5

## test case that does not converge
stderr_dump = redirect_stderr()
@test isnan(estimate_pi0(p, CensoredBUM(0.5, 0.05, 1e-6, 2)))

@test issubtype(typeof(CensoredBUM()), Pi0Estimator)
@test issubtype(typeof(CensoredBUM(0.2, 0.1)), Pi0Estimator)

@test_throws DomainError CensoredBUM(-0.5, 0.05)
@test_throws DomainError CensoredBUM(1.5, 0.05)
@test_throws DomainError CensoredBUM(0.5, -0.05)
@test_throws DomainError CensoredBUM(0.5, 1.05)

@test_throws DomainError CensoredBUM(0.5, 0.05, -1e-6, 100)
@test_throws DomainError CensoredBUM(0.5, 0.05, 1.1, 100)
@test_throws DomainError CensoredBUM(0.5, 0.05, 1e-6, -10)

f = fit(CensoredBUM(), p)
@test issubtype(typeof(f), CensoredBUMFit)
@test_approx_eq_eps f.π0 0.55797 1e-5

pi0est, pars, is_converged = MultipleTesting.cbum_pi0(ones(50))
@test_approx_eq pi0est 1.0
@test is_converged

pi0est, pars, is_converged = MultipleTesting.cbum_pi0_naive(ones(50))
@test isnan(pi0est)
@test !is_converged


## BUM_pi0 ##
## needs better test cases and reference values
println(" ** ", "BUM_pi0")

@test_approx_eq_eps estimate_pi0(p, BUM(0.5)) 0.55528 1e-5
@test_approx_eq_eps estimate_pi0(p, BUM()) 0.55528 1e-5
@test_approx_eq_eps estimate_pi0(p0, BUM()) 1.0 1e-5
@test_approx_eq_eps estimate_pi0(p1, BUM()) 0.10874 1e-5

## test case that does not converge
@test isnan(estimate_pi0(p, BUM(0.5, 1e-6, 2)))

@test issubtype(typeof(BUM()), Pi0Estimator)

@test_throws DomainError BUM(-0.5)
@test_throws DomainError BUM(1.5)


## Flat grenander ##
println(" ** ", "FlatGrenander_pi0")

FlatGrenander = MultipleTesting.FlatGrenander

@test issubtype(typeof(FlatGrenander()), Pi0Estimator)

pu = collect(0.1:0.05:0.9)
@test_approx_eq estimate_pi0(pu, FlatGrenander()) 1.0
@test_approx_eq estimate_pi0(pu.^0.5, FlatGrenander()) 1.0

@test_approx_eq estimate_pi0(p0, FlatGrenander()) 1.0
@test estimate_pi0(p1, FlatGrenander()) < 0.15
@test_approx_eq_eps estimate_pi0(p, FlatGrenander()) pi0 0.1

## longest constant interval: low level
lci = MultipleTesting.longest_constant_interval

p = [0:0.1:1;];
f = [1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.1];
@test_approx_eq lci(p, f) 0.5

f = [1.5, 1.5, 1.5, 1.5, 1.5, 1.2, 1.2, 1.2, 0.2, 0.2, 0.1];
@test_approx_eq lci(p, f) 0.2

f = [0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1];
@test_approx_eq lci(p, f) 0.2

f = [0.5, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1];
@test_approx_eq lci(p, f) 0.5

p = [0.1, 0.3, 0.5, 0.9];
f = [0.5, 0.5, 0.2, 0.2];
@test_approx_eq lci(p, f) 0.2

p = [0.1, 0.5, 0.7, 0.9];
f = [0.5, 0.5, 0.2, 0.2];
@test_approx_eq lci(p, f) 0.5

end
