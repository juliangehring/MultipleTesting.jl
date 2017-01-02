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
    validPValues(pValues)
    k = length(pValues)
    x = -2 * sum(log(pValues))
    p = ccdf(Chisq(2k), x)
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
    validPValues(pValues)
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
    validPValues(pValues)
    k = length(pValues)
    z = cquantile(Normal(), pValues)
    z = sum(z) ./ sqrt(k)
    p = ccdf(Normal(), z)
    return p
end

function stouffer_combination(pValues, weights)
    validPValues(pValues)
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

function tippett_combination(pValues)
    validPValues(pValues)
    k = length(pValues)
    p = 1.0 - (1.0 - minimum(pValues))^k
    return p
end
