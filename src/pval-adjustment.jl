## p-value adjustment methods ##

function bonferroni{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    return min(pValues * length(pValues), 1.)
end


function benjamini_hochberg{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, bejamini_hochberg_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

bejamini_hochberg_multiplier(i::Int, n::Int) = n/(n-i)


function benjamini_hochberg{T<:FloatingPoint}(pValues::Vector{T}, pi0::T)
    validPValues([pi0])
    benjamini_hochberg(pValues) .* pi0
end


function holm{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepdown!(sortedPValues, holm_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

holm_multiplier(i::Int, n::Int) = (n-i+1)


function benjamini_yekutieli{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, benjamini_yekutieli_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

function benjamini_yekutieli_multiplier(i::Int, n::Int)
    c = sum([1/i for i in 1:n])
    return ((n*c)/(n-i))
end


function hochberg{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, hochberg_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

hochberg_multiplier(i::Int, n::Int) = (i+1)


function hommel{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    q = fill(minimum(n .* pValues./[1:n; ]), n)
    pa = fill(q[1], n)
    for j in (n-1):-1:2
        ij = 1:(n-j+1)
        i2 = (n-j+2):n
        q1 = minimum(j .* sortedPValues[i2]./([2:j; ]))
        q[ij] = min(j .* sortedPValues[ij], q1)
        q[i2] = q[n-j+1]
        pa = max(pa, q)
    end
    max(pa, sortedPValues)[originalOrder]
end
