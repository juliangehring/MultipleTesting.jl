## pi_0 estimators ##

## https://github.com/StoreyLab/qvalue/blob/master/R/pi0est.R
function storey_pi0{T<:FloatingPoint}(pValues::Vector{T}, lambda::T)
    pi0 = (sum(pValues .>= lambda)) / (1-lambda) / length(pValues)
    pi0 = min(pi0, 1.)
    return pi0
end


function bootstrap_pi0{T<:FloatingPoint}(pValues::Vector{T}, lambda::Vector{T} = [0.05:0.05:0.95], q::T = 0.1)
    #validPvalues(pValues)
    #validPvalues(lambda) ## TODO check bounds
    n = length(pValues)
    ## CHCK: Only for smoothing (due to cubic spline)?
    if length(lambda) < 4
        throw(ArgumentError())
    end
    lambda = sort(lambda)
    pi0 = [mean(pValues .>= l) / (1-l) for l in lambda]
    min_pi0 = quantile(pi0, q)
    ## in a loop? relevant only for large vectors 'lambda'
    w = [sum(pValues .>= l) for l in lambda]
    mse = (w ./ (n .^ 2 .* (1-lambda) .^ 2 )) .* (1-w/n) + (pi0-min_pi0) .^2
    pi0 = min(pi0[indmin(mse)], 1.)
    pi0
end

#pval = [0.01:0.05:0.91] ## consistent with 'qvalue::pi0est'
#bootstrap_pi0(pval)
