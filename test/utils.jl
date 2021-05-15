# helper functions

import Random: shuffle!

function unsort(x; kws...)
    y = copy(x)
    while issorted(y; kws...)
        shuffle!(y)
    end
    return y
end

function unorder(x; kws...)
    ord = collect(1:length(x))
    while issorted(x[ord]; kws...)
        shuffle!(ord)
    end
    return ord
end