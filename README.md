# MultipleTesting

The `MultipleTesting` package offers common algorithms for p-value adjustment and π0 estimation.

![xkcd p-value guide](pvalues.png)


## Package Status

[![Package Status](http://pkg.julialang.org/badges/MultipleTesting_0.4.svg)](http://pkg.julialang.org/?pkg=MultipleTesting)
[![Linux/Mac Build Status](https://travis-ci.org/julian-gehring/MultipleTesting.jl.svg?branch=master)](https://travis-ci.org/julian-gehring/MultipleTesting.jl)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/1ld0ppptisirryt1/branch/master?svg=true)](https://ci.appveyor.com/project/julian-gehring/multipletesting-jl/branch/master)
[![Coverage Status](http://codecov.io/github/julian-gehring/MultipleTesting.jl/coverage.svg?branch=master)](http://codecov.io/github/julian-gehring/MultipleTesting.jl?branch=master&view=all)
[![Documentation Status](https://readthedocs.org/projects/multipletestingjl/badge/?version=master)](http://multipletestingjl.readthedocs.org/en/docs/)


## Features

### p-values adjustment

* Bonferroni
* Benjamini-Hochberg
* Benjamini-Hochberg "Adaptive" with known π0 or π0 estimation method (see section below)
* Benjamini-Yekutieli
* Benjamini-Liu
* Hochberg
* Holm
* Hommel
* Sidak
* ForwardStop

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
```


### π0 estimation

* Storey
* Storey's closed-form bootstrap
* Least SLope (LSL)
* Two STep (TST)
* RightBoundary (Storey's estimate with dynamically chosen λ)
* Censored BUM
* BUM
* FlatGrenander

```julia
estimate_pi0(pvals, Storey())
estimate_pi0(pvals, StoreyBootstrap())
estimate_pi0(pvals, LeastSlope())
estimate_pi0(pvals, TwoStep())
estimate_pi0(pvals, TwoStep(0.05))
estimate_pi0(pvals, RightBoundary())
estimate_pi0(pvals, CensoredBUM())
estimate_pi0(pvals, BUM())
estimate_pi0(pvals, FlatGrenander())
```
