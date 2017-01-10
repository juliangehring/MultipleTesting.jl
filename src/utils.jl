## utility functions ##

function stepup!{T<:AbstractFloat}(sortedPValues::Vector{T}, multiplier::Function, n::Integer = length(sortedPValues))
    sortedPValues[n] *= multiplier(0, n)
    for i in 1:(n-1)
        sortedPValues[n-i] = min(sortedPValues[n-i+1], sortedPValues[n-i] * multiplier(i, n))
    end
    return sortedPValues
end

# multiplier stepdown
function stepdown!{T<:AbstractFloat}(sortedPValues::Vector{T}, multiplier::Function, n::Integer = length(sortedPValues))
  stepfun(p::T, i::Int, n::Int) = p * multiplier(i, n)
  general_stepdown!(sortedPValues, stepfun, n)
  return sortedPValues
end

function general_stepdown!{T<:AbstractFloat}(sortedPValues::Vector{T}, stepfun::Function, n::Integer = length(sortedPValues))
    sortedPValues[1] = stepfun(sortedPValues[1], 1, n)
    for i in 2:n
        sortedPValues[i] = max(sortedPValues[i-1], stepfun(sortedPValues[i], i, n))
    end
    return sortedPValues
end


function reorder{T<:Number}(values::Vector{T})
    newOrder = sortperm(values)
    oldOrder = sortperm(newOrder)
    return newOrder, oldOrder
end


function validPValues{T<:AbstractFloat}(x::Vector{T})
    if !isin(x)
        throw(DomainError())
    end
end


function isin(x::Real, lower::Real = 0., upper::Real = 1.)
    x >= lower && x <= upper
end

function isin{T<:Real}(x::Vector{T}, lower::Real = 0., upper::Real = 1.)
    ex = extrema(x)
    ex[1] >= lower && ex[2] <= upper
end


function isotonic_regression_reference{T<:AbstractFloat}(y::Vector{T}, w::Vector{T})
    #todo: ignore zero weights
    y = copy(y)
    w = copy(w)
    m = length(y)
    cnts = ones(Int64, m)
    i = 2
    # ... not most efficient way but could be fun to (ab)use iterator protocol
    while !done(y, i)
        if y[i] < y[i-1]
            y[i-1] = (w[i]*y[i]+w[i-1]*y[i-1])/(w[i]+w[i-1])
            w[i-1] = w[i]+w[i-1]
            cnts[i-1] += cnts[i]
            deleteat!(y, i)
            deleteat!(w, i)
            deleteat!(cnts, i)
            i = max(i-2, 1)
        end
        i += 1
    end
    yisotonic = vcat([y[idx]*ones(Float64, cnt) for (idx, cnt) in enumerate(cnts)]...)
    return yisotonic
end

function isotonic_regression_reference{T<:AbstractFloat}(y::Vector{T})
    isotonic_regression_reference(y, ones(y))
end


function isotonic_regression{T<:AbstractFloat}(y::Vector{T}, weights::Vector{T})
    n = length(y)
    if n <= 1
        return y
    end
    if n != length(weights)
        throw(DimensionMismatch("Lengths of values and weights mismatch"))
    end
    @inbounds begin
        n -= 1
        while true
            i = 1
            is_pooled = false
            while i <= n
                k = i
                while k <= n && y[k] >= y[k+1]
                    k += 1
                end
                if y[i] != y[k]
                    numerator = 0.0
                    denominator = 0.0
                    for j in i:k
                        numerator += y[j] * weights[j]
                        denominator += weights[j]
                    end
                    m = numerator / denominator
                    for j in i:k
                        y[j] = m
                    end
                    is_pooled = true
                end
                i = k + 1
            end
            if !is_pooled
               break
            end
        end
    end
    return y
end

function isotonic_regression{T<:AbstractFloat}(y::Vector{T})
    isotonic_regression(y, ones(y))
end


function grenander{T<:AbstractFloat}(pv::Vector{T})
    pv_sorted = sort(pv)
    ## ecdf that handles duplicated values
    pv_sorted_unique, counts = rle(pv_sorted)
    ecdf_value = cumsum(counts)
    ecdf_value = ecdf_value ./ ecdf_value[end]

    Δx = diff(pv_sorted_unique)
    Δy = diff(ecdf_value)

    f = Δy ./ Δx
    f = -isotonic_regression(-f, Δx)
    F  = ecdf_value[1] + vcat(0, cumsum(f .* Δx))
    f = push!(f, f[end])

    return pv_sorted_unique, f, F
end
