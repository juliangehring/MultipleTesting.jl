## qValues ##

function qValues{T<:AbstractFloat}(pValues::AbstractVector{T}, pi0::AbstractFloat, pfdr::Bool = false)
    valid_pvalues(pValues)
    valid_pvalues([pi0])
    n = length(pValues)
    sortedIndex = sortperm(pValues)
    qValues = zeros(pValues)
    q = T(Inf)
    for i in n:-1:1
        idx = sortedIndex[i]
        if pfdr
            q = min(pValues[idx] * n / (i * (1 - (1 - pValues[idx]) ^ n)), q)
        else
            q = min(pValues[idx] * n / i, q)
        end
        qValues[idx] = pi0 * q
    end
    return qValues
end
