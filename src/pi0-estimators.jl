## pi_0 estimators ##

## https://github.com/StoreyLab/qvalue/blob/master/R/pi0est.R
function storey_pi0{T<:FloatingPoint}(pValues::Vector{T}, lambda::T)
    pi0 = (sum(pValues .>= lambda)) / (1-lambda) / length(pValues)
    pi0 = min(pi0, 1.)
    return pi0
end
