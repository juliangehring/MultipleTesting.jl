# MultipleTesting

The `MultipleTesting` package offers common algorithms for p-value adjustment
and combination as well as the estimation of the proportion π₀ of true null
hypotheses.

![xkcd p-value guide](https://imgs.xkcd.com/comics/p_values.png)


## Features

### Adjustment of p-values

* Bonferroni
* Benjamini-Hochberg
* Adaptive Benjamini-Hochberg with known π₀ or π₀ estimation method (see section below)
* Benjamini-Yekutieli
* Benjamini-Liu
* Hochberg
* Holm
* Hommel
* Sidak
* Forward Stop
* Barber-Candès

```julia
adjust(pvals, Bonferroni())
adjust(pvals, BenjaminiHochberg())
adjust(pvals, BenjaminiHochbergAdaptive(0.9))
adjust(pvals, BenjaminiHochbergAdaptive(Storey()))
adjust(pvals, BenjaminiYekutieli())
adjust(pvals, BenjaminiLiu())
adjust(pvals, Hochberg())
adjust(pvals, Holm())
adjust(pvals, Hommel())
adjust(pvals, Sidak())
adjust(pvals, ForwardStop())
adjust(pvals, BarberCandes())
```

The adjustment can also be performed on the `k` smallest out of `n` p-values:

```julia
adjust(pvals, n, PValueAdjustmentMethod)
```


### Estimation of π₀

* Storey
* Storey's closed-form bootstrap
* Least Slope
* Two Step
* RightBoundary (Storey's estimate with dynamically chosen λ)
* Beta-Uniform Mixture (BUM)
* Censored BUM
* Flat Grenander
* Oracle for known π₀

```julia
estimate(pvals, Storey())
estimate(pvals, StoreyBootstrap())
estimate(pvals, LeastSlope())
estimate(pvals, TwoStep())
estimate(pvals, TwoStep(0.05))
estimate(pvals, TwoStep(0.05, BenjaminiHochbergAdaptive(0.9))
estimate(pvals, RightBoundary())
estimate(pvals, CensoredBUM())
estimate(pvals, BUM())
estimate(pvals, FlatGrenander())
estimate(pvals, Oracle(0.9))
```


### Combination of p-values

* Fisher
* Stouffer, optionally with weights
* Logit
* Tippett
* Simes
* Wilkinson
* Minimum of adjusted p-values

```julia
combine(pvals, FisherCombination())
combine(pvals, StoufferCombination())
combine(pvals, weights, StoufferCombination())
combine(pvals, LogitCombination())
combine(pvals, TippettCombination())
combine(pvals, SimesCombination())
combine(pvals, WilkinsonCombination(rank))
combine(pvals, MinimumCombination(PValueAdjustment()))
```


### Higher criticism

* Higher criticism scores
* Higher criticism threshold

```julia
estimate(pvals, HigherCriticismScores())
estimate(pvals, HigherCriticismThreshold())
```


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
