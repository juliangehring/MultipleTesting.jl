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

[![xkcd p-value guide](https://imgs.xkcd.com/comics/p_values.png)](https://xkcd.com/license.html "XKCD")


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
combine(pvals, Fisher())
combine(pvals, Stouffer())
combine(pvals, weights, Stouffer())
combine(pvals, Logit())
combine(pvals, Tippett())
combine(pvals, Simes())
combine(pvals, Wilkinson(rank))
combine(pvals, Minimum(PValueAdjustment()))
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
