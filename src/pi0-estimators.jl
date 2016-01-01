### estimators for π0 (pi0) ###

abstract Pi0Estimator


## Storey estimator ##

type Storey <: Pi0Estimator
    λ::AbstractFloat

    Storey(λ) = isin(λ, 0., 1.) ? new(λ) : throw(DomainError())
end

Storey() = Storey(0.1)

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::Storey)
    storey_pi0(pValues, pi0estimator.λ)
end

function storey_pi0{T<:AbstractFloat}(pValues::Vector{T}, lambda::T)
    validPValues(pValues)
    pi0 = (sum(pValues .>= lambda) / length(pValues)) / (1.-lambda)
    pi0 = min(pi0, 1.)
    return pi0
end


## Storey bootstrap estimator ##

type StoreyBootstrap <: Pi0Estimator
    ## check range of arguments
    λseq::Vector{AbstractFloat}
    q   ::AbstractFloat

    StoreyBootstrap(λseq, q) =
        isin(λseq, 0., 1.) && isin(q, 0., 1.) ? new(λseq, q) : throw(DomainError())
end

StoreyBootstrap() = StoreyBootstrap(collect(0.05:0.05:0.95), 0.1)

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::StoreyBootstrap)
    bootstrap_pi0(pValues, pi0estimator.λseq, pi0estimator.q)
end

function bootstrap_pi0{T<:AbstractFloat,S<:AbstractFloat}(pValues::Vector{T}, lambda::Vector{S} = collect(0.05:0.05:0.95), q::S = 0.1)
    #validPValues(pValues)
    n = length(pValues)
    if !issorted(lambda)
        lambda = sort(lambda)
    end
    pi0 = Float64[mean(pValues .>= l) / (1-l) for l in lambda]
    min_pi0 = quantile(pi0, q)
    ## in a loop? relevant only for very large vectors 'lambda'
    w = Int[sum(pValues .>= l) for l in lambda]
    mse = (w ./ (n .^ 2 .* (1-lambda) .^ 2 )) .* (1-w/n) + (pi0-min_pi0) .^2
    pi0 = min(pi0[indmin(mse)], 1.)
    pi0
end


## Least SLope (LSL) estimator ##

type LeastSlope <: Pi0Estimator
end

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::LeastSlope)
    lsl_pi0(pValues)
end

function lsl_pi0{T<:AbstractFloat}(pValues::Vector{T})
    n = length(pValues)
    if !issorted(pValues)
        pValues = sort(pValues)
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

## alternative, vectorized version
## used for comparison and compactness
function lsl_pi0_vec{T<:AbstractFloat}(pValues::Vector{T})
    n = length(pValues)
    if !issorted(pValues)
        pValues = sort(pValues)
    end
    s = (1 - pValues) ./ (n - collect(1:n) + 1)
    d = diff(s) .< 0
    idx = findfirst(d) + 1
    pi0 = min( 1/s[idx] + 1, n ) / n
    return(pi0)
end


## Oracle

type Oracle <: Pi0Estimator
    π0::AbstractFloat

    Oracle(π0) = isin(π0, 0., 1.) ? new(π0) : throw(DomainError())
end

Oracle() = Oracle(1.0)

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::Oracle)
    pi0estimator.π0
end
