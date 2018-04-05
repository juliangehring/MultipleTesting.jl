## qValues ##

struct QValues
    pi0estimator::Pi0Estimator
    pfdr::Bool
end

QValues() = QValues(Oracle(1.0), false)


function estimate(pValues::PValues{T}, method::QValues) where T<:AbstractFloat
    pi0 = estimate_pi0(pValues, method.pi0estimator)
    qValues(pValues, pi0, method.pfdr)
end

function qValues(pValues::PValues{T}, pi0::AbstractFloat, pfdr::Bool) where T<:AbstractFloat
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
