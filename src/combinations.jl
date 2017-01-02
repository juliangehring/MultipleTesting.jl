### Combination methods for p-values ###


## Fisher combination ##

type FisherCombination <: PValueCombinationMethod
end

function combine{T<:AbstractFloat}(pValues::Vector{T}, method::FisherCombination)
    fisher_combination(pValues)
end

function fisher_combination(pValues)
    validPValues(pValues)
    k = length(pValues)
    if k == 1
        return pValues[1]
    end
    if minimum(pValues) == 0.0
        return NaN
    end
    x = -2 * sum(log(pValues)) # p = 0 => Inf
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
    validPValues(pValues)
    k = length(pValues)
    if k == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0.0 || pmax == 1.0
        return NaN
    end
    c = sqrt( (5k+2)*k*pi^2 / ((5k+4)*3) )
    x = -sum(log(pValues./(1-pValues))) / c
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
    if k == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0.0 || pmax == 1.0
        return NaN
    end
    z = cquantile(Normal(), pValues)
    z = sum(z) ./ sqrt(k)
    p = ccdf(Normal(), z)
    return p
end

function stouffer_combination(pValues, weights)
    validPValues(pValues)
    k = length(pValues)
    if k == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0.0 || pmax == 1.0
        return NaN
    end
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
    if k == 1
        return pValues[1]
    end
    p = 1.0 - (1.0 - minimum(pValues))^k
    return p
end
