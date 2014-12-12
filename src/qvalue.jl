## qValues ##

## https://github.com/StoreyLab/qvalue/blob/master/R/qvalue.R
function qValues{T<:FloatingPoint}(pValues::Vector{T}, pi0::T, pfdr::Bool = false)
    n = length(pValues)
    u = sortperm(pValues)
    v = competerank(pValues) ## ties with 'min'
    ##v = rank(p, ties.method="max") 
    if pfdr
        qvals = (pi0 .* n .* pValues) ./ (v .* (1 - (1 - pValues) .^ n))
    else
        qvals = (pi0 .* n .* pValues) ./ v
    end
    qvals[u[n]] = min(qvals[u[n]], 1)
    for i in (n - 1):-1:1
        qvals[u[i]] = min(qvals[u[i]], qvals[u[i + 1]])
    end
    return qvals
end
