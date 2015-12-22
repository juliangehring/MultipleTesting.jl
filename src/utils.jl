## utility functions ##

function stepup!{T<:AbstractFloat}(sortedPValues::Vector{T}, multiplier::Function, n::Integer = length(sortedPValues))
    sortedPValues[n] *= multiplier(0, n)
    for i in 1:(n-1)
        sortedPValues[n-i] = min(sortedPValues[n-i+1], sortedPValues[n-i] * multiplier(i, n))
    end
    return sortedPValues
end


function stepdown!{T<:AbstractFloat}(sortedPValues::Vector{T}, multiplier::Function, n::Integer = length(sortedPValues))
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


function validPValues{T<:AbstractFloat}(x::Vector{T})
    ex = extrema(x)
    if ex[1] .< 0.0 || ex[2] .> 1.0
        throw(DomainError())
    end
end
