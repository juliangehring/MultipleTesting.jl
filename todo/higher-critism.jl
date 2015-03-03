## HC

function hc_score(pValues)
    validPvalues(pValues)
    n = length(pValues)
    f = competerank(pValues) / n ## ties = max
    var = f .* (1-f) / n
    var[var == 0] = min( var[var > 0] ) # just to make sure we have no zero variance
    hc = abs(f-pValues) / sqrt(var)
    hc
end


