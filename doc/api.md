# MultipleTesting


## Functions [Exported]

---

<a id="function__betauniformmixturemodel.1" class="lexicon_definition"></a>
#### MultipleTesting.BetaUniformMixtureModel
Beta Uniform Mixture Model (BUM)

**Arguments**

- π0 : Contributing fraction of the uniform distribution to the full model
- α, β : Parameters of the Beta distribution, Float64, default: 0.5, 3.0

**Returns**

`MixtureModel`, as defined in the `Distributions` package, composed of

- a uniform distribution in the interval [0, 1], with weight/prior π0
- a Beta distribution with parameters α and β, with weight/prior 1-π0

**Examples**

```julia
bum = BetaUniformMixtureModel(0.2)
bum = BetaUniformMixtureModel(0.2, 0.5, 1.0)

using Distributions

rand(bum, 10)

pdf(bum, 0:0.05:1)

cdf(bum, 0:0.05:1)
```



## Types [Exported]

---

<a id="type__bum.1" class="lexicon_definition"></a>
#### MultipleTesting.BUM
BUM π0 estimator

BUM(γ0, xtol, maxiter)


---

<a id="type__censoredbum.1" class="lexicon_definition"></a>
#### MultipleTesting.CensoredBUM
Censored BUM π0 estimator

CensoredBUM(γ0, λ, xtol, maxiter)


---

<a id="type__leastslope.1" class="lexicon_definition"></a>
#### MultipleTesting.LeastSlope
Least SLope (LSL) π0 estimator

LeastSlope()


---

<a id="type__oracle.1" class="lexicon_definition"></a>
#### MultipleTesting.Oracle
Oracle π0

Oracle(π0)


---

<a id="type__rightboundary.1" class="lexicon_definition"></a>
#### MultipleTesting.RightBoundary
Right boundary π0 estimator

RightBoundary(λseq)


---

<a id="type__storey.1" class="lexicon_definition"></a>
#### MultipleTesting.Storey
Storey π0 estimator

**Parameters**

- λ : tuning parameter, FloatingPoint, default: 0.1

**Examples**

```julia
Storey()
Storey(0.1)
```

**References**

Storey, JD (2002). "A Direct Approach to False Discovery Rates." Journal of the
Royal Statistical Society, doi:10.1111/1467-9868.00346



---

<a id="type__storeybootstrap.1" class="lexicon_definition"></a>
#### MultipleTesting.StoreyBootstrap
Storey closed-form bootstrap π0 estimator

StoreyBootstrap(λseq, q)

Reference: David Robinson, 2012


---

<a id="type__twostep.1" class="lexicon_definition"></a>
#### MultipleTesting.TwoStep
Two step π0 estimator

TwoStep(α)

Reference: Benjamini, Krieger and Yekutieli, 2006


