abstract Pi0Estimator

abstract Pi0Fit

abstract GeneralizedTypeOneError

type FWER <: GeneralizedTypeOneError
end

type FDR <: GeneralizedTypeOneError
end

abstract PValueAdjustmentMethod{typeI} <: GeneralizedTypeOneError

