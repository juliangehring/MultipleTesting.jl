## Higher criticism

immutable HigherCriticismScores end

function estimate{T<:AbstractFloat}(pValues::PValues{T}, method::HigherCriticismScores)
    higher_criticism_scores(pValues)
end

function higher_criticism_scores{T<:AbstractFloat}(pValues::PValues{T})
    n = length(pValues)
    F = (n+1 - competerank(-pValues)) ./ n  # ECDF
    denom = F .* (1.0 - F) ./ n
    # avoid denominator of 0 for last value
    idx0 = denom .== 0
    denom[idx0] = minimum(denom[.!idx0]) + eps()  # conservative
    hcs = abs.(F - pValues) ./ sqrt.(denom)
    return hcs
end


immutable HigherCriticalValue end

function estimate{T<:AbstractFloat}(pValues::PValues{T}, method::HigherCriticalValue)
    higher_critical_value(pValues)
end

function higher_critical_value{T<:AbstractFloat}(pValues::PValues{T})
    idx_hcv = indmax(higher_criticism_scores(pValues))
    return pValues[idx_hcv]
end
