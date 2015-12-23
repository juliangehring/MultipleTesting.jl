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


function isotonicregression(y::Array{Float64,1},w::Array{Float64,1})
  #todo: ignore zero weights
  y=copy(y)
  w=copy(w)
  m = length(y)
  cnts = ones(Int64,m)
  i = 2
  # ... not most efficient way but could be fun to (ab)use iterator protocol
  while (!done(y,i))
    if y[i]<y[i-1]
      y[i-1]=(w[i]*y[i]+w[i-1]*y[i-1])/(w[i]+w[i-1])
      w[i-1]=w[i]+w[i-1]
      cnts[i-1] += cnts[i]
      deleteat!(y,i)
      deleteat!(w,i)
      deleteat!(cnts,i)
      i = max(i-2,1)
    end
    i += 1
  end
  yisotonic = vcat([y[idx]*ones(Float64,cnt) for (idx,cnt) in enumerate(cnts)]...)
end

function isotonicregression(y::Array{Float64,1})
  isotonicregression(y, ones(Float64, length(y)))
end
