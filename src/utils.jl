## utility functions ##

function stepup!(sortedPValues, multiplier)
    for i in 0:(n - 1)
        if i == 0
            sortedPValues[n-i] *= multiplier
        else
            sortedPValues[n-i] = min(sortedPValues[n-i+1], sortedPValues[n-i] * multiplier)
        end
    end
    return sortedPValues
end


function stepdown!(sortedPValues, multiplier)
    for i in 1:n
      if i == 1
        sortedPValues[i] *= multiplier
      else
        sortedPValues[i] = max(sortedPValues[i-1], sortedPValues[i] * multiplier)
      end
    end
    return sortedPValues
end


function reorder(values)
    newOrder = sortperm(values)
    oldOrder = sortperm(newOrder)
    return newOrder, oldOrder
end
