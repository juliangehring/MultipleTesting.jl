# Conservative Estimator
type ConservativeEstimator <: Pi0Estimator
end

const conservative = ConservativeEstimator()

function estimatepi0(pvalues::Array{Float64,1}, pi0estimator::ConservativeEstimator)
  one(Float64)
end
