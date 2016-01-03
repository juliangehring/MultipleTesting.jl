## Test pvalue adjustment methods ##
module Test_pval_adjustment

using MultipleTesting
using Base.Test

pi0 = 0.4
pval1 = [0.0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.75, 1.0]
ref1 = Dict([
             (Bonferroni, [0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0]),
             (Holm, [0.0, 9e-4, 8e-3, 0.07, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0]),
             (Hochberg, [0.0, 9e-4, 8e-3, 0.07, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0]),
             (Hommel, [0.0, 9e-4, 8e-3, 0.07, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0]),
             (BenjaminiHochberg, [0.0, 0.0005, 0.003333333, 0.025, 0.1, 0.166666667, 0.285714286, 0.5, 0.833333333, 1.0]),
             (BenjaminiYekutieli, [0.0, 0.001464484, 0.009763228, 0.073224206, 0.292896825, 0.488161376, 0.836848073, 1.0, 1.0, 1.0]),
             (Sidak, [ 0.0, 0.0009995501, 0.0099551198, 0.0956179250, 0.4012630608, 0.6513215599, 0.8926258176, 0.9939533824, 0.9999990463, 1.0000000000]) #mutoss
             ])

pval2 = [0.0001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.1, 0.4, 0.75, 1.0]
ref2 = Dict([
             (Bonferroni, [0.001, 0.001, 0.01, 0.1, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0]),
             (Holm, [0.001, 0.001, 0.008, 0.07, 0.3, 0.5, 0.5, 1.0, 1.0, 1.0]),
             (Hochberg, [9e-4, 9e-4, 8e-3, 0.07, 0.3, 0.4, 0.4, 1.0, 1.0, 1.0]),
             (Hommel, [9e-4, 9e-4, 8e-3, 0.07, 0.25, 0.4, 0.4, 1.0, 1.0, 1.0]),
             (BenjaminiHochberg, [0.0005, 0.0005, 0.003333333, 0.025, 0.1, 0.142857143, 0.142857143, 0.5, 0.833333333, 1.0]),
             (BenjaminiYekutieli, [0.001464484, 0.001464484, 0.009763228, 0.073224206, 0.292896825, 0.418424036, 0.418424036, 1.0, 1.0, 1.0]),
             (Sidak,  [0.0009995501, 0.0009995501, 0.0099551198, 0.0956179250, 0.4012630608, 0.6513215599, 0.6513215599, 0.9939533824, 0.9999990463, 1.0000000000]) #mutoss
             ])

for method in keys(ref1)
    println(" ** ", method)
    @test issubtype(method, PValueAdjustmentMethod)
    @test issubtype(typeof(method()), PValueAdjustmentMethod)
    @test_throws MethodError method(0.1)
    ## no valid p-values as input 
    @test_throws DomainError adjust([-1.0, 0.5], method())
    @test_throws DomainError adjust([0.5, 1.5], method())
    ## single p-value is returned unchanged
    pval = rand(1)
    @test adjust(pval, method()) == pval
    ## compare with reference values
    @test_approx_eq_eps adjust(pval1, method()) ref1[method] 1e-9
    ## compare with reference values having ties
    @test_approx_eq_eps adjust(pval2, method()) ref2[method] 1e-9
end

end

