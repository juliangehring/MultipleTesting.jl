## p-value adjustment methods ##

function bonferroni{T<:FloatingPoint}(pValues::Vector{T})
    return min(pValues * length(pValues), 1.)
end


function benjamini_hochberg{T<:FloatingPoint}(pValues::Vector{T})
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


function holm{T<:FloatingPoint}(pValues::Vector{T})
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
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    #c = sum([1/n for i in 1:n])
    stepup!(sortedPValues, benjamini_yekutieli_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

function bejamini_yekutieli_multiplier(i::Int, n::Int)
    c = sum([1/i for i in 1:n])
    return ((n*c)/(n-i))
end


function hochberg{T<:FloatingPoint}(pValues::Vector{T})
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
