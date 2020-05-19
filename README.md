# MultipleTesting

The `MultipleTesting` package offers common algorithms for p-value adjustment
and combination as well as the estimation of the proportion π₀ of true null
hypotheses.

![xkcd p-value guide](https://imgs.xkcd.com/comics/p_values.png)


## Features


### Adjustment of p-values

```julia
adjust(pvalues, <:PValueAdjustmentMethod)
```

The adjustment can also be performed on the `k` smallest out of `n` p-values:

```julia
adjust(pvalues, n, <:PValueAdjustmentMethod)
```


#### Bonferroni
    
```julia
adjust(pvalues, Bonferroni())
```

Bonferroni, C.E. (1936). Teoria statistica delle classi e calcolo delle probabilita
(Libreria internazionale Seeber).


#### Benjamini-Hochberg

```julia
adjust(pvalues, BenjaminiHochberg())
```

Adaptive Benjamini-Hochberg with known π₀ or π₀ estimation method (see section below)

```julia
adjust(pvalues, BenjaminiHochbergAdaptive(π₀))
adjust(pvalues, BenjaminiHochbergAdaptive(<:PValueAdjustmentMethod))
```

Benjamini, Y., and Hochberg, Y. (1995). Controlling the False Discovery Rate: A
Practical and Powerful Approach to Multiple Testing. Journal of the Royal
Statistical Society. Series B (Methodological) 57, 289–300.


#### Benjamini-Yekutieli

```julia
adjust(pvalues, BenjaminiYekutieli())
```

Benjamini, Y., and Yekutieli, D. (2001). The Control of the False Discovery Rate
in Multiple Testing under Dependency. The Annals of Statistics 29, 1165–1188.


#### Benjamini-Liu

```julia
adjust(pvalues, BenjaminiLiu())
```

Benjamini, Y., and Liu, W. (1999). A step-down multiple hypotheses testing
procedure that controls the false discovery rate under independence. Journal of
Statistical Planning and Inference 82, 163–170.


#### Hochberg

```julia
adjust(pvalues, Hochberg())
```

Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of
significance. Biometrika 75, 800–802.


#### Holm

```julia
adjust(pvalues, Holm())
```

Holm, S. (1979). A Simple Sequentially Rejective Multiple Test Procedure.
Scandinavian Journal of Statistics 6, 65–70.


#### Hommel

```julia
adjust(pvalues, Hommel())
```

Hommel, G. (1988). A stagewise rejective multiple test procedure based on a
modified Bonferroni test. Biometrika 75, 383–386.


#### Sidak

```julia
adjust(pvalues, Sidak())
```

Šidák, Z. (1967). Rectangular Confidence Regions for the Means of Multivariate
Normal Distributions. Journal of the American Statistical Association 62,
626–633.


#### Forward Stop

```julia
adjust(pvalues, ForwardStop())
```

G’Sell, M.G., Wager, S., Chouldechova, A., and Tibshirani, R. (2016). Sequential
selection procedures and false discovery rate control. J. R. Stat. Soc. B 78,
423–444.


#### Barber-Candès

```julia
adjust(pvalues, BarberCandes())
```

Barber, R.F., and Candès, E.J. (2015). Controlling the false discovery rate via
knockoffs. Ann. Statist. 43, 2055–2085.

Arias-Castro, E., and Chen, S. (2017). Distribution-free multiple testing.
Electron. J. Statist. 11, 1983–2001.



### Estimation of π₀

```julia
estimate(pvalues, <:Pi0Estimator)
```


#### Storey

```julia
estimate(pvalues, Storey())
```

Storey, J.D., Taylor, J.E., and Siegmund, D. (2004). Strong control,
conservative point estimation and simultaneous conservative consistency of false
discovery rates: a unified approach. Journal of the Royal Statistical Society:
Series B (Statistical Methodology) 66, 187–205.


#### Storey's closed-form bootstrap

```julia
estimate(pvalues, StoreyBootstrap())
```

Robinson, D. (2016). Original Procedure for Choosing λ.
http://varianceexplained.org/files/pi0boot.pdf


#### Least Slope

```julia
estimate(pvalues, LeastSlope())
```

Benjamini, Y., and Hochberg, Y. (2000). On the Adaptive Control of the False
Discovery Rate in Multiple Testing With Independent Statistics. Journal of
Educational and Behavioral Statistics 25, 60–83.


#### Two Step

```julia
estimate(pvalues, TwoStep())
estimate(pvalues, TwoStep(α))
estimate(pvalues, TwoStep(α, <:PValueAdjustmentMethod)
```

Benjamini, Y., Krieger, A.M., and Yekutieli, D. (2006). Adaptive linear step-up
procedures that control the false discovery rate. Biometrika 93, 491–507.


#### RightBoundary

Storey's estimate with dynamically chosen λ

```julia
estimate(pvalues, RightBoundary())
```

Liang, K., and Nettleton, D. (2012). Adaptive and dynamic adaptive procedures
for false discovery rate control and estimation. Journal of the Royal
Statistical Society: Series B (Statistical Methodology) 74, 163–182.


#### Beta-Uniform Mixture (BUM)

```julia
estimate(pvalues, BUM())
```

Pounds, S., and Morris, S.W. (2003). Estimating the occurrence of false
positives and false negatives in microarray studies by approximating and
partitioning the empirical distribution of p-values. Bioinformatics 19,
1236–1242.


#### Censored BUM

```julia
estimate(pvalues, CensoredBUM())
```

Markitsis, A., and Lai, Y. (2010). A censored beta mixture model for the
estimation of the proportion of non-differentially expressed genes.
Bioinformatics 26, 640–646.


#### Flat Grenander

```julia
estimate(pvalues, FlatGrenander())
```

Langaas, M., Lindqvist, B.H., and Ferkingstad, E. (2005). Estimating the
proportion of true null hypotheses, with application to DNA microarray data.
Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67,
555–572.


#### ConvexDecreasing

```julia
estimate(pvalues, ConvexDecreasing())
fit(ConvexDecreasing(), pvalues)
```

Langaas, M., Lindqvist, B.H., and Ferkingstad, E. (2005). Estimating the
proportion of true null hypotheses, with application to DNA microarray data.
Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67,
555–572.


#### Oracle for known π₀

```julia
estimate(pvalues, Oracle(π₀))
```



### Combination of p-values

```julia
combine(pvalues, <:PValueCombination)
```


#### Fisher

```julia
combine(pvalues, FisherCombination())
```

Fisher, R.A. (1925). Statistical methods for research workers (Genesis
Publishing Pvt Ltd).


#### Stouffer

Optionally with weights

```julia
combine(pvalues, StoufferCombination())
combine(pvalues, weights, StoufferCombination())
```

Stouffer, S.A. (1949). The American soldier. Vol. 1: Adjustment during army life
(Princeton University Press).

Liptak, T. (1958). On the combination of independent tests. Magyar Tud Akad Mat
Kutato Int Kozl 3, 171–197.


#### Logit

```julia
combine(pvalues, LogitCombination())
```

Mudholkar, G.S., and George, E.O. (1977). The Logit Statistic for Combining
Probabilities - An Overview (Rochester University NY, Dept of Statistics).


#### Tippett

```julia
combine(pvalues, TippettCombination())
```

Tippett, L.H.C. (1931). The Methods of Statistics. An introduction mainly for
workers in the biological sciences.


#### Simes

```julia
combine(pvalues, SimesCombination())
```

Simes, R.J. (1986). An improved Bonferroni procedure for multiple tests of
significance. Biometrika 73, 751–754.


#### Wilkinson

```julia
combine(pvalues, WilkinsonCombination(rank))
```

Wilkinson, B. (1951). A statistical consideration in psychological research.
Psychological Bulletin 48, 156.


#### Minimum of adjusted p-values

```julia
combine(pvalues, MinimumCombination(PValueAdjustment()))
```



### Higher criticism

Higher criticism scores and threshold

```julia
estimate(pvalues, HigherCriticismScores())
estimate(pvalues, HigherCriticismThreshold())
```

Donoho, D., and Jin, J. (2008). Higher criticism thresholding: Optimal feature
selection when useful features are rare and weak. PNAS 105, 14790–14795.

Klaus, B., and Strimmer, K. (2013). Signal identification for rare and weak
features: higher criticism or false discovery rates? Biostatistics 14, 129–143.



### Simulation of p-value distributions


#### Beta Uniform Mixture Model

```julia
BetaUniformMixtureModel(π₀, α, β)
```

Pounds, S., and Morris, S.W. (2003). Estimating the occurrence of false
positives and false negatives in microarray studies by approximating and
partitioning the empirical distribution of p-values. Bioinformatics 19,
1236–1242.


## Installation

The `MultipleTesting` package is part of the Julia ecosphere and the latest
release version can be installed with

```julia
pkg> add MultipleTesting
```

More details on packages and how to manage them can be found in the
[package section](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html)
of the Julia documentation.


## Feedback and Contributions

Contributions of any kind are very welcome. Please feel free to open pull
requests or issues with your questions or ideas.


## Package Status

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliangehring.github.io/MultipleTesting.jl/stable)

[![DOI](https://zenodo.org/badge/27935122.svg)](https://zenodo.org/badge/latestdoi/27935122)

![Testing](https://github.com/juliangehring/MultipleTesting.jl/workflows/Testing/badge.svg)

[![Coverage](https://codecov.io/gh/juliangehring/MultipleTesting.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliangehring/MultipleTesting.jl)

The package uses [semantic versioning](https://semver.org/).
