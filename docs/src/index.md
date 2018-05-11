# MultipleTesting.jl

```@meta
CurrentModule = MultipleTesting
DocTestSetup = quote
    using MultipleTesting
end
```


## Package Features

## Manual

```@contents
```

## Library

### Adjustment of p-Values

```@docs
Bonferroni
BenjaminiHochberg
BenjaminiHochbergAdaptive
BenjaminiYekutieli
BenjaminiLiu
Hochberg
Holm
Hommel
Sidak
ForwardStop
BarberCandes
```

### Estimation of π₀

```@docs
Storey
StoreyBootstrap
LeastSlope
Oracle
TwoStep
RightBoundary
CensoredBUM
BUM
FlatGrenander
ConvexDecreasing
```

### Combination of p-Values

```@docs
FisherCombination
LogitCombination
StoufferCombination
TippettCombination
SimesCombination
WilkinsonCombination
MinimumCombination
```

### Higher Criticism

```@docs
HigherCriticismScores
HigherCriticismThreshold
```

### Modelling of p-Value Distributions

```@docs
BetaUniformMixtureModel
```


## Index

```@index
```
