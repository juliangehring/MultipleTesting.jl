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
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    if minimum(pValues) == 0
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
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0 || pmax == 1
        return NaN
    end
    c = sqrt( (5n+2)*n*pi^2 / ((5n+4)*3) )
    x = -sum(log.(pValues./(1 .- pValues))) / c
    p = ccdf(TDist(5n+4), x)
    return p
end


## Stouffer combination ##

struct StoufferCombination <: PValueCombination
end

function combine(pValues::PValues{T}, method::StoufferCombination) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0 || pmax == 1
        return NaN
    end
    z = cquantile.(Normal(), pValues)
    z = sum(z) ./ sqrt(n)
    p = ccdf(Normal(), z)
    return p
end

function combine(pValues::PValues{T}, weights::Vector{T}, method::StoufferCombination) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0 || pmax == 1
        return NaN
    end
    z = cquantile.(Normal(), pValues) .* weights
    z = sum(z) ./ sqrt(sum(abs2, weights))
    p = ccdf(Normal(), z)
    return p
end

function combine(pValues::PValues{T}, weights::Weights, method::StoufferCombination) where T<:AbstractFloat
    combine(pValues, values(weights), method)
end


## Tippett combination ##

struct TippettCombination <: PValueCombination
end

function combine(pValues::PValues{T}, method::TippettCombination) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    p = 1 - (1 - minimum(pValues))^n
    return p
end


## Simes combination ##

struct SimesCombination <: PValueCombination
end

function combine(pValues::PValues{T}, method::SimesCombination) where T<:AbstractFloat
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
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    rank = method.rank
    if rank < 1 || rank > n
        throw(ArgumentError("Rank must be in 1,..,$(n)"))
    end
    p_rank = sort(pValues)[rank]
    p = cdf(Beta(rank, n-rank+1), p_rank)
    return p
end


## Generalised minimum combination ##

struct MinimumCombination <: PValueCombination
    adjustment::PValueAdjustment
end

function combine(pValues::PValues{T}, method::MinimumCombination) where T<:AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    padj = adjust(pValues, method.adjustment)
    p = minimum(padj)
    return p
end
