## Test pvalue adjustment methods ##

using MultipleTesting
using Base.Test

methods = ( bonferroni, holm, hochberg, benjamini_hochberg )

for m in methods
    @test_throws MethodError m()
    ## no valid p-values as input 
    #@test_throws DomainError m(randn(100))
    ## single p-value is returned unchanged
    pval = rand(1)
    @test m(pval) == pval
end
