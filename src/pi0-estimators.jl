## pi_0 estimators ##

abstract Pi0Estimator



# Conservative Estimator
type ConservativeEstimator <: Pi0Estimator
end

const conservative = ConservativeEstimator()

function estimatepi0(pvalues::Array{Float64,1}, pi0estimator::ConservativeEstimator)
  one(Float64)
end

# Storey Estimator
type StoreyEstimator <: Pi0Estimator
  λ::Float64 # actually in [0,1]
end

function estimatepi0(pvalues::Array{Float64,1}, pi0estimator::StoreyEstimator)
  validPValues(pvalues)
  λ = pi0estimator.λ
  m = length(pvalues)
  min(1,sum(pvalues .>= λ)/(1-λ)/length(pvalues))
end

## https://github.com/StoreyLab/qvalue/blob/master/R/pi0est.R
function storey_pi0{T<:AbstractFloat}(pValues::Vector{T}, lambda::T)
    validPValues(pValues)
    pi0 = (sum(pValues .>= lambda) / length(pValues)) / (1-lambda)
    pi0 = min(pi0, 1.)
    return pi0
end

# Storey Bootstrap Estimator

type StoreyBootstrapEstimator <: Pi0Estimator
  λseq::Vector{Float64}
  q   ::Float64
end

StoreyBootstrapEstimator() = StoreyBootstrapEstimator([0.05:0.05:0.95], 0.1)

function estimatepi0(pvalues::Array{Float64,1}, pi0estimator::StoreyBootstrapEstimator)
  bootstrap_pi0(pvalues, pi0estimator.λseq, pi0estimator.q)
end


function bootstrap_pi0{T<:AbstractFloat}(pValues::Vector{T}, lambda::Vector{T} = collect(0.05:0.05:0.95), q::T = 0.1)
    #validPValues(pValues)
    #validPValues(lambda) ## TODO check bounds
    n = length(pValues)
    if !issorted(lambda)
        sort!(lambda)
    end
    pi0 = Float64[mean(pValues .>= l) / (1-l) for l in lambda]
    min_pi0 = quantile(pi0, q)
    ## in a loop? relevant only for very large vectors 'lambda'
    w = Int[sum(pValues .>= l) for l in lambda]
    mse = (w ./ (n .^ 2 .* (1-lambda) .^ 2 )) .* (1-w/n) + (pi0-min_pi0) .^2
    pi0 = min(pi0[indmin(mse)], 1.)
    pi0
end


function lsl_pi0_vec{T<:AbstractFloat}(pValues::Vector{T})
    n = length(pValues)
    ## sorting requires most time
    if !issorted(pValues)
        sort!(pValues)
    end
    s = (1 - pValues) ./ (n - collect(1:n) + 1)
    d = diff(s) .< 0
    idx = findfirst(d) + 1
    pi0 = min( 1/s[idx] + 1, n ) / n
    return(pi0)
end


function lsl_pi0{T<:AbstractFloat}(pValues::Vector{T})
    n = length(pValues)
    ## sorting requires most time
    if !issorted(pValues)
        sort!(pValues)
    end
    s0 = lsl_slope(1, n, pValues)
    sx = 0.
    for i in 2:n
        s1 = lsl_slope(i, n, pValues)
        if (s1 - s0) < 0.
            sx = s1
            break
        end
        s0 = s1
    end
    pi0 = min( 1/sx + 1, n ) / n
    return(pi0)
end

function lsl_slope{T<:AbstractFloat}(i::Int, n::Int, pval::Vector{T})
    s = (1 - pval[i]) / (n - i + 1)
    return s
end
