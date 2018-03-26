### Combination methods for p-values ###

function combine(pValues::AbstractVector{T}, method::M) where {T<:AbstractFloat, M<:PValueCombination}
    combine(PValues(pValues), method)
end

function combine(pValues::AbstractVector{T}, weights::Weights, method::M) where {T<:AbstractFloat, M<:PValueCombination}
    combine(PValues(pValues), weights, method)
end

function combine(pValues::AbstractVector{T}, weights::AbstractVector{R}, method::M) where {T<:AbstractFloat, R<:Real, M<:PValueCombination}
    combine(PValues(pValues), weights, method)
end


## Fisher combination ##

struct FisherCombination <: PValueCombination
end

function combine(pValues::PValues{T}, method::FisherCombination) where T<:AbstractFloat
    fisher_combination(pValues)
end

function fisher_combination(pValues::PValues{T}) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    if minimum(pValues) == 0.0
        return NaN
    end
    x = -2 * sum(log.(pValues))
    p = ccdf(Chisq(2n), x)
    return p
end


## Logit combination ##

struct LogitCombination <: PValueCombination
end

function combine(pValues::PValues{T}, method::LogitCombination) where T<:AbstractFloat
    logit_combination(pValues)
end

function logit_combination(pValues::PValues{T}) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0.0 || pmax == 1.0
        return NaN
    end
    c = sqrt( (5n+2)*n*pi^2 / ((5n+4)*3) )
    x = -sum(log.(pValues./(1-pValues))) / c
    p = ccdf(TDist(5n+4), x)
    return p
end


## Stouffer combination ##

struct StoufferCombination <: PValueCombination
end

function combine(pValues::PValues{T}, method::StoufferCombination) where T<:AbstractFloat
    stouffer_combination(pValues)
end

function combine(pValues::PValues{T}, weights::Vector{T}, method::StoufferCombination) where T<:AbstractFloat
    stouffer_combination(pValues, weights)
end

function combine(pValues::PValues{T}, weights::Weights, method::StoufferCombination) where T<:AbstractFloat
    stouffer_combination(pValues, values(weights))
end

function stouffer_combination(pValues::PValues{T}) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0.0 || pmax == 1.0
        return NaN
    end
    z = cquantile.(Normal(), pValues)
    z = sum(z) ./ sqrt(n)
    p = ccdf(Normal(), z)
    return p
end

function stouffer_combination(pValues::PValues{T}, weights::Vector{T}) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0.0 || pmax == 1.0
        return NaN
    end
    z = cquantile.(Normal(), pValues) .* weights
    z = sum(z) ./ sqrt(sum(weights.^2))
    p = ccdf(Normal(), z)
    return p
end


## Tippett combination ##

struct TippettCombination <: PValueCombination
end

function combine(pValues::PValues{T}, method::TippettCombination) where T<:AbstractFloat
    tippett_combination(pValues)
end

function tippett_combination(pValues::PValues{T}) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    p = 1.0 - (1.0 - minimum(pValues))^n
    return p
end


## Simes combination ##

struct SimesCombination <: PValueCombination
end

function combine(pValues::PValues{T}, method::SimesCombination) where T<:AbstractFloat
    simes_combination(pValues)
end

function simes_combination(pValues::PValues{T}) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pValues = sort(pValues)  # faster than `sortperm`
    p = n * minimum(pValues./(1:n))
    return p
end


## Wilkinson combination ##

struct WilkinsonCombination <: PValueCombination
    rank::Integer

    function WilkinsonCombination(rank)
        if rank < 1
            throw(ArgumentError("Rank must be positive."))
        end
        return new(rank)
    end
end

function combine(pValues::PValues{T}, method::WilkinsonCombination) where T<:AbstractFloat
    wilkinson_combination(pValues, method.rank)
end

function wilkinson_combination(pValues::PValues{T}, rank::Integer) where T<:AbstractFloat
    n = length(pValues)
    if rank < 1 || rank > n
        throw(ArgumentError("Rank must be in 1,..,$(n)"))
    end
    if n == 1
        return pValues[1]
    end
    p_rank = sort(pValues)[rank]
    p = cdf(Beta(rank, n-rank+1), p_rank)
    return p
end


## Generalised minimum combination ##

struct MinimumCombination <: PValueCombination
    method::PValueAdjustment
end

function combine(pValues::PValues{T}, method::MinimumCombination) where T<:AbstractFloat
    minimum_combination(pValues, method.method)
end

function minimum_combination(pValues::PValues{T}, pAdjustMethod::PValueAdjustment) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    padj = adjust(pValues, pAdjustMethod)
    p = minimum(padj)
    return p
end
