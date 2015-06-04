## Test pvalue adjustment methods ##
module Test_pval_adjustment

using MultipleTesting
using Base.Test

pval1 = [0.0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.75, 1.0]
ref1 = Dict([
             (bonferroni, [0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0]),
             (holm, [0.0, 9e-4, 8e-3, 0.07, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0]),
             (hochberg, [0.0, 9e-4, 8e-3, 0.07, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0]),
             (hommel, [0.0, 9e-4, 8e-3, 0.07, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0]),
             (benjamini_hochberg, [0.0, 0.0005, 0.003333333, 0.025, 0.1, 0.166666667, 0.285714286, 0.5, 0.833333333, 1.0]),
             (benjamini_yekutieli, [0.0, 0.001464484, 0.009763228, 0.073224206, 0.292896825, 0.488161376, 0.836848073, 1.0, 1.0, 1.0])
             ])

pval2 = [0.0001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.1, 0.4, 0.75, 1.0]
ref2 = Dict([
             (bonferroni, [0.001, 0.001, 0.01, 0.1, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0]),
             (holm, [0.001, 0.001, 0.008, 0.07, 0.3, 0.5, 0.5, 1.0, 1.0, 1.0]),
             (hochberg, [9e-4, 9e-4, 8e-3, 0.07, 0.3, 0.4, 0.4, 1.0, 1.0, 1.0]),
             (hommel, [9e-4, 9e-4, 8e-3, 0.07, 0.25, 0.4, 0.4, 1.0, 1.0, 1.0]),
             (benjamini_hochberg, [0.0005, 0.0005, 0.003333333, 0.025, 0.1, 0.142857143, 0.142857143, 0.5, 0.833333333, 1.0]),
             (benjamini_yekutieli, [0.001464484, 0.001464484, 0.009763228, 0.073224206, 0.292896825, 0.418424036, 0.418424036, 1.0, 1.0, 1.0])
             ])

for m in keys(ref1)
    println(" ** ", m)
    @test_throws MethodError m()
    ## no integers as input
    @test_throws MethodError m([0, 1])
    ## no valid p-values as input 
    @test_throws DomainError m([-1.0, 0.5])
    @test_throws DomainError m([0.5, 1.5])
    ## single p-value is returned unchanged
    pval = rand(1)
    @test m(pval) == pval
    ## compare with reference values
    @test_approx_eq_eps m(pval1) ref1[m] 1e-9
    ## compare with reference values having ties
    @test_approx_eq_eps m(pval2) ref2[m] 1e-9
end

end
