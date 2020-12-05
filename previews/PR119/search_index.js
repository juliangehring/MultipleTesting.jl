var documenterSearchIndex = {"docs":
[{"location":"combination/#Combination-of-p-Values-1","page":"Combination of p-Values","title":"Combination of p-Values","text":"","category":"section"},{"location":"combination/#","page":"Combination of p-Values","title":"Combination of p-Values","text":"CurrentModule = MultipleTesting\nDocTestSetup = quote\n    using MultipleTesting\nend","category":"page"},{"location":"combination/#Methods-1","page":"Combination of p-Values","title":"Methods","text":"","category":"section"},{"location":"combination/#","page":"Combination of p-Values","title":"Combination of p-Values","text":"combine","category":"page"},{"location":"combination/#MultipleTesting.combine","page":"Combination of p-Values","title":"MultipleTesting.combine","text":"combine(PValues, <:PValueCombination)\n\nCombine p-values\n\nExamples\n\njulia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);\n\njulia> combine(pvals, FisherCombination())\n0.007616871850449092\njulia> combine(pvals, StoufferCombination())\n0.007098326181265917\n\nSee also\n\nPValueCombinations:\n\nFisherCombination LogitCombination StoufferCombination TippettCombination SimesCombination WilkinsonCombination MinimumCombination\n\n\n\n\n\n","category":"function"},{"location":"combination/#Types-1","page":"Combination of p-Values","title":"Types","text":"","category":"section"},{"location":"combination/#","page":"Combination of p-Values","title":"Combination of p-Values","text":"FisherCombination\nLogitCombination\nStoufferCombination\nTippettCombination\nSimesCombination\nWilkinsonCombination\nMinimumCombination","category":"page"},{"location":"combination/#MultipleTesting.FisherCombination","page":"Combination of p-Values","title":"MultipleTesting.FisherCombination","text":"Fisher's p-value combination\n\nExamples\n\njulia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);\n\njulia> combine(pvals, FisherCombination())\n0.007616871850449092\n\nReferences\n\nFisher, R.A. (1925). Statistical methods for research workers (Genesis Publishing Pvt Ltd).\n\n\n\n\n\n","category":"type"},{"location":"combination/#MultipleTesting.LogitCombination","page":"Combination of p-Values","title":"MultipleTesting.LogitCombination","text":"Logit p-value combination\n\nExamples\n\njulia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);\n\njulia> combine(pvals, LogitCombination())\n0.006434494635148462\n\nReferences\n\nMudholkar, G.S., and George, E.O. (1977). The Logit Statistic for Combining Probabilities - An Overview (Rochester University NY, Dept of Statistics).\n\n\n\n\n\n","category":"type"},{"location":"combination/#MultipleTesting.StoufferCombination","page":"Combination of p-Values","title":"MultipleTesting.StoufferCombination","text":"Stouffer's p-value combination\n\nExamples\n\njulia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);\n\njulia> combine(pvals, StoufferCombination())\n0.007098326181265917\n\njulia> weights = [1.0, 2.0, 0.4, 1.5];\n\njulia> combine(pvals, weights, StoufferCombination())\n0.007331653763696742\n\nReferences\n\nStouffer, S.A. (1949). The American soldier. Vol. 1: Adjustment during army life (Princeton University Press).\n\nLiptak, T. (1958). On the combination of independent tests. Magyar Tud Akad Mat Kutato Int Kozl 3, 171–197.\n\n\n\n\n\n","category":"type"},{"location":"combination/#MultipleTesting.TippettCombination","page":"Combination of p-Values","title":"MultipleTesting.TippettCombination","text":"Tippett's p-value combination\n\nExamples\n\njulia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);\n\njulia> combine(pvals, TippettCombination())\n0.039403990000000055\n\nReferences\n\nTippett, L.H.C. (1931). The Methods of Statistics. An introduction mainly for workers in the biological sciences.\n\n\n\n\n\n","category":"type"},{"location":"combination/#MultipleTesting.SimesCombination","page":"Combination of p-Values","title":"MultipleTesting.SimesCombination","text":"Simes's p-value combination\n\nExamples\n\njulia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);\n\njulia> combine(pvals, SimesCombination())\n0.04\n\nReferences\n\nSimes, R.J. (1986). An improved Bonferroni procedure for multiple tests of significance. Biometrika 73, 751–754.\n\n\n\n\n\n","category":"type"},{"location":"combination/#MultipleTesting.WilkinsonCombination","page":"Combination of p-Values","title":"MultipleTesting.WilkinsonCombination","text":"Wilkinson's p-value combination\n\nExamples\n\njulia> pv = PValues([0.01, 0.02, 0.3, 0.5]);\n\njulia> combine(pv, WilkinsonCombination(1))  # combination with rank 1\n0.03940399000000003\n\njulia> combine(pv, WilkinsonCombination(4))  # combination with rank 4\n0.0625\n\nReferences\n\nWilkinson, B. (1951). A statistical consideration in psychological research. Psychological Bulletin 48, 156.\n\n\n\n\n\n","category":"type"},{"location":"combination/#MultipleTesting.MinimumCombination","page":"Combination of p-Values","title":"MultipleTesting.MinimumCombination","text":"Minimum of adjusted p-value combination\n\nExamples\n\njulia> pv = PValues([0.01, 0.02, 0.3, 0.5]);\n\njulia> combine(pv, MinimumCombination(BenjaminiHochberg()))\n0.04\n\njulia> combine(pv, MinimumCombination(ForwardStop()))\n0.01005033585350145\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#Adjustment-of-p-Values-1","page":"Adjustment of p-Values","title":"Adjustment of p-Values","text":"","category":"section"},{"location":"adjustment/#","page":"Adjustment of p-Values","title":"Adjustment of p-Values","text":"CurrentModule = MultipleTesting\nDocTestSetup = quote\n    using MultipleTesting\nend","category":"page"},{"location":"adjustment/#Methods-1","page":"Adjustment of p-Values","title":"Methods","text":"","category":"section"},{"location":"adjustment/#","page":"Adjustment of p-Values","title":"Adjustment of p-Values","text":"adjust","category":"page"},{"location":"adjustment/#MultipleTesting.adjust","page":"Adjustment of p-Values","title":"MultipleTesting.adjust","text":"adjust(PValues, <:PValueAdjustment)\nadjust(PValues, Int, <:PValueAdjustment)\n\nAdjustment of p-values\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, BenjaminiHochberg())\n4-element Array{Float64,1}:\n 0.004\n 0.02\n 0.039999999999999994\n 0.5\n\njulia> adjust(pvals, 6, BenjaminiHochberg()) # 4 out of 6 p-values\n4-element Array{Float64,1}:\n 0.006\n 0.03\n 0.06\n 0.75\n\njulia> adjust(pvals, BarberCandes())\n4-element Array{Float64,1}:\n 0.3333333333333333\n 0.3333333333333333\n 0.3333333333333333\n 1.0\n\nSee also\n\nPValueAdjustments:\n\nBonferroni BenjaminiHochberg BenjaminiHochbergAdaptive BenjaminiYekutieli BenjaminiLiu Hochberg Holm Hommel Sidak ForwardStop BarberCandes\n\n\n\n\n\n","category":"function"},{"location":"adjustment/#Types-1","page":"Adjustment of p-Values","title":"Types","text":"","category":"section"},{"location":"adjustment/#","page":"Adjustment of p-Values","title":"Adjustment of p-Values","text":"Bonferroni\nBenjaminiHochberg\nBenjaminiHochbergAdaptive\nBenjaminiYekutieli\nBenjaminiLiu\nHochberg\nHolm\nHommel\nSidak\nForwardStop\nBarberCandes","category":"page"},{"location":"adjustment/#MultipleTesting.Bonferroni","page":"Adjustment of p-Values","title":"MultipleTesting.Bonferroni","text":"Bonferroni adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, Bonferroni())\n4-element Array{Float64,1}:\n 0.004\n 0.04\n 0.12\n 1.0\n\njulia> adjust(pvals, 6, Bonferroni())\n4-element Array{Float64,1}:\n 0.006\n 0.06\n 0.18\n 1.0\n\nReferences\n\nBonferroni, C.E. (1936). Teoria statistica delle classi e calcolo delle probabilita (Libreria internazionale Seeber).\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#MultipleTesting.BenjaminiHochberg","page":"Adjustment of p-Values","title":"MultipleTesting.BenjaminiHochberg","text":"Benjamini-Hochberg adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, BenjaminiHochberg())\n4-element Array{Float64,1}:\n 0.004\n 0.02\n 0.039999999999999994\n 0.5\n\njulia> adjust(pvals, 6, BenjaminiHochberg())\n4-element Array{Float64,1}:\n 0.006\n 0.03\n 0.06\n 0.75\n\nReferences\n\nBenjamini, Y., and Hochberg, Y. (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society. Series B (Methodological) 57, 289–300.\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#MultipleTesting.BenjaminiHochbergAdaptive","page":"Adjustment of p-Values","title":"MultipleTesting.BenjaminiHochbergAdaptive","text":"Adaptive Benjamini-Hochberg adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, BenjaminiHochbergAdaptive(Oracle(0.5))) # known π₀ of 0.5\n4-element Array{Float64,1}:\n 0.002\n 0.01\n 0.019999999999999997\n 0.25\n\njulia> adjust(pvals, BenjaminiHochbergAdaptive(StoreyBootstrap())) # π₀ estimator\n4-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n\njulia> adjust(pvals, 6, BenjaminiHochbergAdaptive(StoreyBootstrap()))\n4-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n\nReferences\n\nBenjamini, Y., Krieger, A. M. & Yekutieli, D. (2006). Adaptive linear step-up  procedures that control the false discovery rate. Biometrika 93, 491–507.\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#MultipleTesting.BenjaminiYekutieli","page":"Adjustment of p-Values","title":"MultipleTesting.BenjaminiYekutieli","text":"Benjamini-Yekutieli adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, BenjaminiYekutieli())\n4-element Array{Float64,1}:\n 0.00833333333333334\n 0.0416666666666667\n 0.0833333333333334\n 1.0\n\njulia> adjust(pvals, 6, BenjaminiYekutieli())\n4-element Array{Float64,1}:\n 0.01470000000000001\n 0.07350000000000005\n 0.14700000000000008\n 1.0\n\nReferences\n\nBenjamini, Y., and Yekutieli, D. (2001). The Control of the False Discovery Rate in Multiple Testing under Dependency. The Annals of Statistics 29, 1165–1188.\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#MultipleTesting.BenjaminiLiu","page":"Adjustment of p-Values","title":"MultipleTesting.BenjaminiLiu","text":"Benjamini-Liu adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, BenjaminiLiu())\n4-element Array{Float64,1}:\n 0.003994003998999962\n 0.022275749999999983\n 0.02955000000000002\n 0.125\n\njulia> adjust(pvals, 6, BenjaminiLiu())\n4-element Array{Float64,1}:\n 0.0059850199850060015\n 0.04084162508333341\n 0.07647146000000005\n 0.4375\n\nReferences\n\nBenjamini, Y., and Liu, W. (1999). A step-down multiple hypotheses testing procedure that controls the false discovery rate under independence. Journal of Statistical Planning and Inference 82, 163–170.\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#MultipleTesting.Hochberg","page":"Adjustment of p-Values","title":"MultipleTesting.Hochberg","text":"Hochberg adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, Hochberg())\n4-element Array{Float64,1}:\n 0.004\n 0.03\n 0.06\n 0.5\n\njulia> adjust(pvals, 6, Hochberg())\n4-element Array{Float64,1}:\n 0.006\n 0.05\n 0.12\n 1.0\n\nReferences\n\nHochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika 75, 800–802.\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#MultipleTesting.Holm","page":"Adjustment of p-Values","title":"MultipleTesting.Holm","text":"Holm adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, Holm())\n4-element Array{Float64,1}:\n 0.004\n 0.03\n 0.06\n 0.5\n\njulia> adjust(pvals, 6, Holm())\n4-element Array{Float64,1}:\n 0.006\n 0.05\n 0.12\n 1.0\n\nReferences\n\nHolm, S. (1979). A Simple Sequentially Rejective Multiple Test Procedure. Scandinavian Journal of Statistics 6, 65–70.\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#MultipleTesting.Hommel","page":"Adjustment of p-Values","title":"MultipleTesting.Hommel","text":"Hommel adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, Hommel())\n4-element Array{Float64,1}:\n 0.004\n 0.03\n 0.06\n 0.5\n\njulia> adjust(pvals, 6, Hommel())\n4-element Array{Float64,1}:\n 0.006\n 0.05\n 0.12\n 1.0\n\nReferences\n\nHommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika 75, 383–386.\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#MultipleTesting.Sidak","page":"Adjustment of p-Values","title":"MultipleTesting.Sidak","text":"Šidák adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, Sidak())\n4-element Array{Float64,1}:\n 0.003994003998999962\n 0.039403990000000055\n 0.11470719000000007\n 0.9375\n\njulia> adjust(pvals, 6, Sidak())\n4-element Array{Float64,1}:\n 0.0059850199850060015\n 0.058519850599\n 0.1670279950710002\n 0.984375\n\nReferences\n\nŠidák, Z. (1967). Rectangular Confidence Regions for the Means of Multivariate Normal Distributions. Journal of the American Statistical Association 62, 626–633.\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#MultipleTesting.ForwardStop","page":"Adjustment of p-Values","title":"MultipleTesting.ForwardStop","text":"Forward-Stop adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, ForwardStop())\n4-element Array{Float64,1}:\n 0.0010005003335835344\n 0.005525418093542492\n 0.013836681223931188\n 0.1836643060579347\n\njulia> adjust(pvals, 6, ForwardStop())\n4-element Array{Float64,1}:\n 0.0010005003335835344\n 0.005525418093542492\n 0.013836681223931188\n 0.1836643060579347\n\nReferences\n\nG’Sell, M.G., Wager, S., Chouldechova, A., and Tibshirani, R. (2016). Sequential selection procedures and false discovery rate control. J. R. Stat. Soc. B 78, 423–444.\n\n\n\n\n\n","category":"type"},{"location":"adjustment/#MultipleTesting.BarberCandes","page":"Adjustment of p-Values","title":"MultipleTesting.BarberCandes","text":"Barber-Candès adjustment\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> adjust(pvals, BarberCandes())\n4-element Array{Float64,1}:\n 0.3333333333333333\n 0.3333333333333333\n 0.3333333333333333\n 1.0\n\nReferences\n\nBarber, R.F., and Candès, E.J. (2015). Controlling the false discovery rate via knockoffs. Ann. Statist. 43, 2055–2085.\n\nArias-Castro, E., and Chen, S. (2017). Distribution-free multiple testing. Electron. J. Statist. 11, 1983–2001.\n\n\n\n\n\n","category":"type"},{"location":"pi0/#Estimation-of-π-1","page":"Estimation of π₀","title":"Estimation of π₀","text":"","category":"section"},{"location":"pi0/#","page":"Estimation of π₀","title":"Estimation of π₀","text":"CurrentModule = MultipleTesting\nDocTestSetup = quote\n    using MultipleTesting\nend","category":"page"},{"location":"pi0/#Methods-1","page":"Estimation of π₀","title":"Methods","text":"","category":"section"},{"location":"pi0/#","page":"Estimation of π₀","title":"Estimation of π₀","text":"estimate","category":"page"},{"location":"pi0/#Distributions.estimate","page":"Estimation of π₀","title":"Distributions.estimate","text":"estimate(PValues, <:Pi0Estimator)\n\nEstimate π₀, the fraction of tests under the null hypothesis\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, StoreyBootstrap())\n0.0\njulia> estimate(pvals, FlatGrenander())\n0.42553191489361697\n\nSee also\n\nPi0Estimators:\n\nStorey StoreyBootstrap LeastSlope Oracle TwoStep RightBoundary CensoredBUM BUM FlatGrenander ConvexDecreasing\n\n\n\n\n\n","category":"function"},{"location":"pi0/#Types-1","page":"Estimation of π₀","title":"Types","text":"","category":"section"},{"location":"pi0/#","page":"Estimation of π₀","title":"Estimation of π₀","text":"Storey\nStoreyBootstrap\nLeastSlope\nOracle\nTwoStep\nRightBoundary\nCensoredBUM\nBUM\nFlatGrenander\nConvexDecreasing","category":"page"},{"location":"pi0/#MultipleTesting.Storey","page":"Estimation of π₀","title":"MultipleTesting.Storey","text":"Storey's π₀ estimator\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, Storey())\n0.22222222222222224\n\njulia> estimate(pvals, Storey(0.4))\n0.33333333333333337\n\nReferences\n\nStorey, J.D., Taylor, J.E., and Siegmund, D. (2004). Strong control, conservative point estimation and simultaneous conservative consistency of false discovery rates: a unified approach. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 66, 187–205.\n\n\n\n\n\n","category":"type"},{"location":"pi0/#MultipleTesting.StoreyBootstrap","page":"Estimation of π₀","title":"MultipleTesting.StoreyBootstrap","text":"Storey's closed-form bootstrap π₀ estimator\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, StoreyBootstrap())\n0.0\n\njulia> estimate(pvals, StoreyBootstrap(0.1:0.1:0.9, 0.2))\n0.0\n\nReferences\n\nRobinson, D. (2016). Original Procedure for Choosing λ. http://varianceexplained.org/files/pi0boot.pdf\n\nStorey, J.D., Taylor, J.E., and Siegmund, D. (2004). Strong control, conservative point estimation and simultaneous conservative consistency of false discovery rates: a unified approach. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 66, 187–205.\n\n\n\n\n\n","category":"type"},{"location":"pi0/#MultipleTesting.LeastSlope","page":"Estimation of π₀","title":"MultipleTesting.LeastSlope","text":"Least Slope (LSL) π₀ estimator\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, LeastSlope())\n1.0\n\nReferences\n\nBenjamini, Y., and Hochberg, Y. (2000). On the Adaptive Control of the False Discovery Rate in Multiple Testing With Independent Statistics. Journal of Educational and Behavioral Statistics 25, 60–83.\n\n\n\n\n\n","category":"type"},{"location":"pi0/#MultipleTesting.Oracle","page":"Estimation of π₀","title":"MultipleTesting.Oracle","text":"Oracle π₀\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, Oracle(0.5)) # a bit boring...\n0.5\n\n\n\n\n\n","category":"type"},{"location":"pi0/#MultipleTesting.TwoStep","page":"Estimation of π₀","title":"MultipleTesting.TwoStep","text":"Two-step π₀ estimator\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, TwoStep())\n0.2\n\njulia> estimate(pvals, TwoStep(0.05, BenjaminiLiu()))\n0.2\n\n\nReferences\n\nBenjamini, Y., Krieger, A.M., and Yekutieli, D. (2006). Adaptive linear step-up procedures that control the false discovery rate. Biometrika 93, 491–507.\n\n\n\n\n\n","category":"type"},{"location":"pi0/#MultipleTesting.RightBoundary","page":"Estimation of π₀","title":"MultipleTesting.RightBoundary","text":"Right boundary π₀ estimator\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, RightBoundary())\n0.2127659574468085\n\njulia> estimate(pvals, RightBoundary(0.1:0.1:0.9))\n0.25\n\nReferences\n\nLiang, K., and Nettleton, D. (2012). Adaptive and dynamic adaptive procedures for false discovery rate control and estimation. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 74, 163–182.\n\n\n\n\n\n","category":"type"},{"location":"pi0/#MultipleTesting.CensoredBUM","page":"Estimation of π₀","title":"MultipleTesting.CensoredBUM","text":"Censored Beta-Uniform Mixture (censored BUM) π₀ estimator\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, CensoredBUM())\n0.21052495526400936\n\n\nReferences\n\nMarkitsis, A., and Lai, Y. (2010). A censored beta mixture model for the estimation of the proportion of non-differentially expressed genes. Bioinformatics 26, 640–646.\n\n\n\n\n\n","category":"type"},{"location":"pi0/#MultipleTesting.BUM","page":"Estimation of π₀","title":"MultipleTesting.BUM","text":"Beta-Uniform Mixture (BUM) π₀ estimator\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, BUM())\n0.22802795505154264\n\nReferences\n\nPounds, S., and Morris, S.W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics 19, 1236–1242.\n\n\n\n\n\n","category":"type"},{"location":"pi0/#MultipleTesting.FlatGrenander","page":"Estimation of π₀","title":"MultipleTesting.FlatGrenander","text":"Flat Grenander π₀ estimator\n\nEstimates π₀ by finding the longest constant interval in the Grenander estimator.\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, FlatGrenander())\n0.42553191489361697\n\nReferences\n\nLangaas, M., Lindqvist, B.H., and Ferkingstad, E. (2005). Estimating the proportion of true null hypotheses, with application to DNA microarray data. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67, 555–572.\n\n\n\n\n\n","category":"type"},{"location":"pi0/#MultipleTesting.ConvexDecreasing","page":"Estimation of π₀","title":"MultipleTesting.ConvexDecreasing","text":"Convex Decreasing π₀ estimator\n\nExamples\n\njulia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, ConvexDecreasing())\n0.013007051336745304\n\nReferences\n\nLangaas, M., Lindqvist, B.H., and Ferkingstad, E. (2005). Estimating the proportion of true null hypotheses, with application to DNA microarray data. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67, 555–572.\n\n\n\n\n\n","category":"type"},{"location":"models/#Modelling-of-p-Value-Distributions-1","page":"Modelling of p-Value Distributions","title":"Modelling of p-Value Distributions","text":"","category":"section"},{"location":"models/#","page":"Modelling of p-Value Distributions","title":"Modelling of p-Value Distributions","text":"CurrentModule = MultipleTesting\nDocTestSetup = quote\n    using MultipleTesting\nend","category":"page"},{"location":"models/#","page":"Modelling of p-Value Distributions","title":"Modelling of p-Value Distributions","text":"BetaUniformMixtureModel","category":"page"},{"location":"models/#MultipleTesting.BetaUniformMixtureModel","page":"Modelling of p-Value Distributions","title":"MultipleTesting.BetaUniformMixtureModel","text":"Beta Uniform Mixture (BUM) Model\n\nArguments\n\nπ0 : Contributing fraction of the uniform distribution to the full model\nα, β : Parameters of the Beta distribution, Float64, default: 0.5, 3.0\n\nReturn values\n\nMixtureModel, as defined in the Distributions package, composed of\n\na uniform distribution in the interval [0, 1], with weight/prior π₀\na Beta distribution with parameters α and β, with weight/prior 1-π₀\n\nExamples\n\njulia> bum = BetaUniformMixtureModel(0.2, 0.5, 1.0);\n\njulia> using Distributions\n\njulia> pdf.(bum, 0.2:0.2:1.0)\n5-element Array{Float64,1}:\n 1.094427190999916\n 0.832455532033676\n 0.7163977794943224\n 0.647213595499958\n 0.6000000000000001\n\n\nReferences\n\nPounds, S., and Morris, S.W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics 19, 1236–1242.\n\n\n\n\n\n","category":"function"},{"location":"#MultipleTesting.jl-1","page":"Home","title":"MultipleTesting.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"CurrentModule = MultipleTesting\nDocTestSetup = quote\n    using MultipleTesting\nend","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The MultipleTesting package offers common algorithms for p-value adjustment and combination as well as the estimation of the proportion π₀ of true null hypotheses.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: xkcd p-value guide)","category":"page"},{"location":"#Features-1","page":"Home","title":"Features","text":"","category":"section"},{"location":"#Adjustment-of-p-values-1","page":"Home","title":"Adjustment of p-values","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Bonferroni\nBenjamini-Hochberg\nAdaptive Benjamini-Hochberg with known π₀ or π₀ estimation method (see section below)\nBenjamini-Yekutieli\nBenjamini-Liu\nHochberg\nHolm\nHommel\nSidak\nForward Stop\nBarber-Candès","category":"page"},{"location":"#","page":"Home","title":"Home","text":"adjust(pvals, Bonferroni())\nadjust(pvals, BenjaminiHochberg())\nadjust(pvals, BenjaminiHochbergAdaptive(0.9))\nadjust(pvals, BenjaminiHochbergAdaptive(Storey()))\nadjust(pvals, BenjaminiYekutieli())\nadjust(pvals, BenjaminiLiu())\nadjust(pvals, Hochberg())\nadjust(pvals, Holm())\nadjust(pvals, Hommel())\nadjust(pvals, Sidak())\nadjust(pvals, ForwardStop())\nadjust(pvals, BarberCandes())","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The adjustment can also be performed on the k smallest out of n p-values:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"adjust(pvals, n, PValueAdjustmentMethod)","category":"page"},{"location":"#Estimation-of-π-1","page":"Home","title":"Estimation of π₀","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Storey\nStorey's closed-form bootstrap\nLeast Slope\nTwo Step\nRightBoundary (Storey's estimate with dynamically chosen λ)\nBeta-Uniform Mixture (BUM)\nCensored BUM\nFlat Grenander\nOracle for known π₀","category":"page"},{"location":"#","page":"Home","title":"Home","text":"estimate(pvals, Storey())\nestimate(pvals, StoreyBootstrap())\nestimate(pvals, LeastSlope())\nestimate(pvals, TwoStep())\nestimate(pvals, TwoStep(0.05))\nestimate(pvals, TwoStep(0.05, BenjaminiHochbergAdaptive(0.9))\nestimate(pvals, RightBoundary())\nestimate(pvals, CensoredBUM())\nestimate(pvals, BUM())\nestimate(pvals, FlatGrenander())\nestimate(pvals, Oracle(0.9))","category":"page"},{"location":"#Combination-of-p-values-1","page":"Home","title":"Combination of p-values","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Fisher\nStouffer, optionally with weights\nLogit\nTippett\nSimes\nWilkinson\nMinimum of adjusted p-values","category":"page"},{"location":"#","page":"Home","title":"Home","text":"combine(pvals, FisherCombination())\ncombine(pvals, StoufferCombination())\ncombine(pvals, weights, StoufferCombination())\ncombine(pvals, LogitCombination())\ncombine(pvals, TippettCombination())\ncombine(pvals, SimesCombination())\ncombine(pvals, WilkinsonCombination(rank))\ncombine(pvals, MinimumCombination(PValueAdjustment()))","category":"page"},{"location":"#Higher-criticism-1","page":"Home","title":"Higher criticism","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Higher criticism scores\nHigher criticism threshold","category":"page"},{"location":"#","page":"Home","title":"Home","text":"estimate(pvals, HigherCriticismScores())\nestimate(pvals, HigherCriticismThreshold())","category":"page"},{"location":"#Manual-1","page":"Home","title":"Manual","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"higher-criticism/#Higher-Criticism-1","page":"Higher Criticism","title":"Higher Criticism","text":"","category":"section"},{"location":"higher-criticism/#","page":"Higher Criticism","title":"Higher Criticism","text":"CurrentModule = MultipleTesting\nDocTestSetup = quote\n    using MultipleTesting\nend","category":"page"},{"location":"higher-criticism/#","page":"Higher Criticism","title":"Higher Criticism","text":"HigherCriticismScores\nHigherCriticismThreshold","category":"page"},{"location":"higher-criticism/#MultipleTesting.HigherCriticismScores","page":"Higher Criticism","title":"MultipleTesting.HigherCriticismScores","text":"Higher criticism scores\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, HigherCriticismScores())\n4-element Array{Float64,1}:\n 1.1500817362257345\n 1.96\n 3.3255375505322444\n 2.3094010767584976\n\n\nReferences\n\nDonoho, D., and Jin, J. (2008). Higher criticism thresholding: Optimal feature selection when useful features are rare and weak. PNAS 105, 14790–14795.\n\nKlaus, B., and Strimmer, K. (2013). Signal identification for rare and weak features: higher criticism or false discovery rates? Biostatistics 14, 129–143.\n\n\n\n\n\n","category":"type"},{"location":"higher-criticism/#MultipleTesting.HigherCriticismThreshold","page":"Higher Criticism","title":"MultipleTesting.HigherCriticismThreshold","text":"Higher criticism threshold\n\nExamples\n\njulia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);\n\njulia> estimate(pvals, HigherCriticismThreshold())\n0.03\n\nReferences\n\nDonoho, D., and Jin, J. (2008). Higher criticism thresholding: Optimal feature selection when useful features are rare and weak. PNAS 105, 14790–14795.\n\nKlaus, B., and Strimmer, K. (2013). Signal identification for rare and weak features: higher criticism or false discovery rates? Biostatistics 14, 129–143.\n\n\n\n\n\n","category":"type"}]
}
