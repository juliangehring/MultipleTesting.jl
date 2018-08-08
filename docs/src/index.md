# MultipleTesting.jl

```@meta
CurrentModule = MultipleTesting
DocTestSetup = quote
    using MultipleTesting
end
```


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
estimate_pi0(pvals, Storey())
estimate_pi0(pvals, StoreyBootstrap())
estimate_pi0(pvals, LeastSlope())
estimate_pi0(pvals, TwoStep())
estimate_pi0(pvals, TwoStep(0.05))
estimate_pi0(pvals, TwoStep(0.05, BenjaminiHochbergAdaptive(0.9))
estimate_pi0(pvals, RightBoundary())
estimate_pi0(pvals, CensoredBUM())
estimate_pi0(pvals, BUM())
estimate_pi0(pvals, FlatGrenander())
estimate_pi0(pvals, Oracle(0.9))
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


## Manual

```@contents
```
