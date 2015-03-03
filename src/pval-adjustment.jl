## p-value adjustment methods ##

function bonferroni{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    return min(pValues * length(pValues), 1.)
end

@doc """
# Benjamini-Hochberg p-value adjustement


## Usage

    pAdjusted = benjamini_hochberg{T<:FloatingPoint}(pValues::Vector{T})

Input arguments:

- `pValues`: Vector of p-values that should be adjusted

Return values:

- `pAdjusted`: Vector of adjusted p-values, matching the input `pValues`.


## References

Benjamini, Y. and Hochberg, Y. (1995):
Controlling the false discovery rate: A practical and powerful approach to multiple testing.
Journal of the Royal Statistical Society

http://en.wikipedia.org/wiki/False_discovery_rate#Benjamini.E2.80.93Hochberg_procedure

## Examples

    pOld = rand(100)
    pNew = benjamini_hochberg(pOld)

""" ->
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
    q = fill(minimum(n * pValues./[1:n]), n)
    pa = fill(q[1], n)
    for j in (n-1):2
        ij = [1:(n-j+1)]
        i2 = [(n-j+2):n]
        q1 = minimum(j * sortedPValues[i2]./([2:j]))
        q[ij] = min(j * sortedPValues[ij], q1)
        q[i2] = q[n-j+1]
        pa = max(pa, q)
    end
    max(pa, sortedPValues)[originalOrder]
end
