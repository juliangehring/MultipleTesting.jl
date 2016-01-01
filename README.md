# MultipleTesting

The `MultipleTesting` package offers common algorithms for p-value adjustment and π0 estimation.

![xkcd p-value guide](pvalues.png)


## Package Status

[![Linux/Mac Build Status](https://travis-ci.org/julian-gehring/MultipleTesting.jl.svg?branch=master)](https://travis-ci.org/julian-gehring/MultipleTesting.jl)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/1ld0ppptisirryt1/branch/master?svg=true)](https://ci.appveyor.com/project/julian-gehring/multipletesting-jl/branch/master)
[![Coverage Status](http://codecov.io/github/julian-gehring/MultipleTesting.jl/coverage.svg?branch=master)](http://codecov.io/github/julian-gehring/MultipleTesting.jl?branch=master&view=all)


## Features

### p-values adjustment

* Bonferroni
* Benjamini-Hochberg
* Benjamini-Hochberg "Oracle" with known π0
* Benjamini-Yekutieli
* Hochberg
* Holm
* Hommel

```julia
adjust(pvals, Bonferroni())
adjust(pvals, BenjaminiHochberg())
adjust(pvals, BenjaminiHochbergOracle(0.9))
adjust(pvals, BenjaminiYekutieli())
adjust(pvals, Hochberg())
adjust(pvals, Holm())
adjust(pvals, Hommel())
```

### π0 estimation

* Storey
* Storey's closed-form bootstrap
* Least SLope (LSL)
* Two Step

```julia
estimate_pi0(pvals, Storey())
estimate_pi0(pvals, StoreyBootstrap())
estimate_pi0(pvals, LeastSlope())
estimate_pi0(pvals, TwoStep(0.05))
```
