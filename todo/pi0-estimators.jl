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


# Storey Bootstrap Estimator

type StoreyBootstrapEstimator <: Pi0Estimator
  λseq::Vector{Float64}
  q   ::Float64
end

StoreyBootstrapEstimator() = StoreyBootstrapEstimator([0.05:0.05:0.95], 0.1)

function estimatepi0(pvalues::Array{Float64,1}, pi0estimator::StoreyBootstrapEstimator)
  bootstrap_pi0(pvalues, pi0estimator.λseq, pi0estimator.q)
end
