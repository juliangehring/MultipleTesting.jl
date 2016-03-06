"""
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

"""
function BetaUniformMixtureModel end

function BetaUniformMixtureModel(π0::Float64, α::Float64 = 0.5, β::Float64 = 1.0)
    if !isin(π0, 0., 1.)
        throw(DomainError())
    end
    MixtureModel([Beta(α, β), Uniform()], [1-π0, π0])
end
