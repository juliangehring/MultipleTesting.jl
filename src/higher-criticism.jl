## Higher criticism

function higher_criticism_scores(pValues::PValues)
    n = length(pValues)
    F = (n+1 - competerank(-pValues)) ./ n  # ECDF
    denom = F .* (1.0 - F) ./ n
    # avoid denominator of 0 for last value
    idx0 = denom .== 0
    denom[idx0] = minimum(denom[!idx0]) + eps()  # conservative
    hcs = abs(F - pValues) ./ sqrt(denom)
    return hcs
end


function higher_critical_value(pValues::PValues)
    idx_hcv = indmax(higher_criticism_scores(pValues))
    return pValues[idx_hcv]
end
