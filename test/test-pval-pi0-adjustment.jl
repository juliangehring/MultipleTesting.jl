## Test pvalue adjustment methods ##
module Test_pval_pi0_adjustment

using MultipleTesting
using Base.Test

pval1 = [0.0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.75, 1.0]
pi0 = 0.4
ref1 = Dict([
             (benjamini_hochberg, [0.0, 0.0005, 0.003333333, 0.025, 0.1, 0.166666667, 0.285714286, 0.5, 0.833333333, 1.0] .* pi0),
             (qValues, [0.0, 0.0002, 0.00133333, 0.01, 0.04, 0.0666667, 0.114286, 0.2, 0.333333, 0.4])
             ])

for m in keys(ref1)
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
    @test_approx_eq_eps m(pval1, pi0) ref1[m] 1e-6
    ## missing: qValue pfdr option, ties
end

end
