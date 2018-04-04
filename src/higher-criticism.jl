## Higher criticism

struct HigherCriticismScores end

function estimate(pValues::PValues{T}, method::HigherCriticismScores) where T<:AbstractFloat
    n = length(pValues)
    F = (n+1 .- competerank(-pValues)) ./ n  # ECDF
    denom = F .* (one(T) .- F) ./ n
    # avoid denominator of 0 for last value
    idx0 = denom .== 0
    denom[idx0] = minimum(denom[.!idx0]) .+ eps()  # conservative
    hcs = abs.(F .- pValues) ./ sqrt.(denom)
    return hcs
end


struct HigherCriticismThreshold end

function estimate(pValues::PValues{T}, method::HigherCriticismThreshold) where T<:AbstractFloat
    idx_hcv = indmax(estimate(pValues, HigherCriticismScores()))
    return pValues[idx_hcv]
end
