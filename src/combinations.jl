### Combination methods for p-values ###

# TODO require at least two p-values as input?

## Fisher combination ##

type FisherCombination <: PValueCombinationMethod
end

function combine{T<:AbstractFloat}(pValues::Vector{T}, method::FisherCombination)
    fisher_combination(pValues)
end

function fisher_combination(pValues)
    # TODO limit to p > 0
    #validPValues(pValues)
    k = length(pValues)
    x = -2 * sum(log(pValues))
    p = ccdf(Chisq(2k), x)
    return p
end

#using Distributions
#p = [0.01, 0.02, 0.1]
#w = [2, 2, 2]
#x = [quantile(Chisq(w[i]), 1-p[i]) for i in 1:length(p)]
#y = sum(x)
#q = ccdf(Chisq(sum(w)), y)
#v = 2.*sum(w).^2./(2.*sum(w))

## Negative Fisher combination ##

type NegativeFisherCombination <: PValueCombinationMethod
end

function combine{T<:AbstractFloat}(pValues::Vector{T}, method::NegativeFisherCombination)
    negative_fisher_combination(pValues)
end

function negative_fisher_combination(pValues)
    p = 1.0 - fisher_combination(1.0-pValues)
    return p
end


## Logit combination ##

type LogitCombination <: PValueCombinationMethod
end

function combine{T<:AbstractFloat}(pValues::Vector{T}, method::LogitCombination)
    logit_combination(pValues)
end

function logit_combination(pValues)
    # TODO limit to p > 0 and p < 1
    #validPValues(pValues)
    k = length(pValues)
    c = sqrt( (5k+2)*k*pi^2 / ((5k+4)*3) )
    x = -sum(log(pValues./(1-pValues))) / c # or name 't'
    p = ccdf(TDist(5k+4), x)
    return p
end


## Stouffer combination ##

type StoufferCombination <: PValueCombinationMethod
end

function combine{T<:AbstractFloat}(pValues::Vector{T}, method::StoufferCombination)
    stouffer_combination(pValues)
end

function combine{T<:AbstractFloat}(pValues::Vector{T}, weights::Vector{T}, method::StoufferCombination)
    stouffer_combination(pValues, weights)
end

function stouffer_combination(pValues)
    k = length(pValues)
    z = cquantile(Normal(), pValues)
    z = sum(z) ./ sqrt(k)
    p = ccdf(Normal(), z)
    return p
end

function stouffer_combination(pValues, weights)
    k = length(pValues)
    z = cquantile(Normal(), pValues) .* weights
    z = sum(z) ./ sqrt(sum(weights.^2))
    p = ccdf(Normal(), z)
    return p
end


## Tippett combination ##

type TippettCombination <: PValueCombinationMethod
end

function combine{T<:AbstractFloat}(pValues::Vector{T}, method::TippettCombination)
    tippett_combination(pValues)
end

"""
    tippett_combination(pValues)

Combine p-values to a global p-value using Tippett's method.

# Examples

```jldoctests
tippett_combination([0.01, 0.2])
# output
0.01990000000000003
```

# Related

- [`TippettCombination`](@ref)

"""
function tippett_combination(pValues)
    k = length(pValues)
    p = 1.0 - (1.0 - minimum(pValues))^k
    return p
end
