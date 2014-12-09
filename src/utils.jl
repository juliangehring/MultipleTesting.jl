## utility functions ##

function stepup!{T<:FloatingPoint}(sortedPValues::Vector{T}, multiplier::Function, n::Integer = length(sortedPValues))
    sortedPValues[n] *= multiplier(0, n)
    for i in 1:(n-1)
        sortedPValues[n-i] = min(sortedPValues[n-i+1], sortedPValues[n-i] * multiplier(i, n))
    end
    return sortedPValues
end


function stepdown!{T<:FloatingPoint}(sortedPValues::Vector{T}, multiplier::Function, n::Integer = length(sortedPValues))
    sortedPValues[1] *= multiplier(1, n)
    for i in 2:n
        sortedPValues[i] = max(sortedPValues[i-1], sortedPValues[i] * multiplier(i, n))
    end
    return sortedPValues
end


function reorder{T<:Number}(values::Vector{T})
    newOrder = sortperm(values)
    oldOrder = sortperm(newOrder)
    return newOrder, oldOrder
end


function validPValues{T<:FloatingPoint}(x::Vector{T})
    ex = extrema(x)
    if ex[1] .< 0.0 || ex[2] .> 1.0
        throw DomainError()
    end
    return res
end

function checkPValues{T<:FloatingPoint}(pValues::Vector{T})
    if any((pValues .< 0.0) | (pValues .> 1.0))
        throw DomainError()
    end
end
