## p-value adjustment methods ##

adjust{T<:AbstractFloat, M<:PValueAdjustmentMethod}(pvals::AbstractVector{T}, method::M) = adjust(PValues(pvals), method)


immutable Bonferroni <: PValueAdjustmentMethod
end

adjust(pvals::PValues, method::Bonferroni) = bonferroni(pvals)


function bonferroni(pValues::PValues)
    return min(pValues * length(pValues), 1.)
end

immutable BenjaminiHochberg <: PValueAdjustmentMethod
end

adjust(pvals::PValues, method::BenjaminiHochberg) = benjamini_hochberg(pvals)

function benjamini_hochberg(pValues::PValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, bejamini_hochberg_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

bejamini_hochberg_multiplier(i::Int, n::Int) = n/(n-i)

immutable BenjaminiHochbergAdaptive <: PValueAdjustmentMethod
    pi0estimator::Pi0Estimator
end

## default to BenjaminiHochberg
BenjaminiHochbergAdaptive{T<:AbstractFloat}(π0::T) = BenjaminiHochbergAdaptive(Oracle(π0))

BenjaminiHochbergAdaptive() = BenjaminiHochbergAdaptive(1.0)

function adjust(pvals::PValues, method::BenjaminiHochbergAdaptive)
  π0 = estimate_pi0(pvals, method.pi0estimator)
  benjamini_hochberg(pvals, π0)
end

function benjamini_hochberg{T<:AbstractFloat}(pValues::PValues, pi0::T)
    validPValues([pi0])
    benjamini_hochberg(pValues) .* pi0
end


immutable BenjaminiYekutieli <: PValueAdjustmentMethod
end

adjust(pvals::PValues, method::BenjaminiYekutieli) = benjamini_yekutieli(pvals)

function benjamini_yekutieli(pValues::PValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, benjamini_yekutieli_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

function benjamini_yekutieli_multiplier(i::Int, n::Int)
    c = sum([1/i for i in 1:n])
    return ((n*c)/(n-i))
end

immutable BenjaminiLiu <: PValueAdjustmentMethod
end

adjust(pvals::PValues, method::BenjaminiLiu) = benjamini_liu(pvals)

function benjamini_liu(pValues::PValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    general_stepdown!(sortedPValues, benjamini_liu_step, n)
    min(sortedPValues[originalOrder], 1.)
end

function benjamini_liu_step{T<:AbstractFloat}(p::T, i::Int, n::Int)
    # a bit more involved because cutoffs at significance α have the form:
    # P_(i) <= 1- [1 - min(1, m/(m-i+1)α)]^{1/(m-i+1)}
    adjp = (1-(1-p)^(n-i+1))*(n-i+1)/n
    if adjp*n/(n-i+1) > 1
        return 0.0
    else
        return adjp
    end
end

immutable Hochberg <: PValueAdjustmentMethod
end

adjust(pvals::PValues, method::Hochberg) = hochberg(pvals)

function hochberg(pValues::PValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, hochberg_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

hochberg_multiplier(i::Int, n::Int) = (i+1)


immutable Holm <: PValueAdjustmentMethod
end

adjust(pvals::PValues, method::Holm) = holm(pvals)

function holm(pValues::PValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepdown!(sortedPValues, holm_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

holm_multiplier(i::Int, n::Int) = (n-i+1)


immutable Hommel <: PValueAdjustmentMethod
end

adjust(pvals::PValues, method::Hommel) = hommel(pvals)

function hommel(pValues::PValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    q = fill(minimum(n .* pValues./[1:n; ]), n)
    pa = fill(q[1], n)
    for j in (n-1):-1:2
        ij = 1:(n-j+1)
        i2 = (n-j+2):n
        q1 = minimum(j .* sortedPValues[i2]./([2:j; ]))
        q[ij] = min(j .* sortedPValues[ij], q1)
        q[i2] = q[n-j+1]
        pa = max(pa, q)
    end
    max(pa, sortedPValues)[originalOrder]
end


immutable Sidak <: PValueAdjustmentMethod
end

adjust(pvals::PValues, method::Sidak) = sidak(pvals)

function sidak(pValues::PValues)
    return min(1-(1-pValues).^length(pValues), 1.)
end

immutable ForwardStop <: PValueAdjustmentMethod
end

adjust(pvals::PValues, method::ForwardStop) = forwardstop(pvals)

function forwardstop(pvalues::PValues)
    n = length(pvalues)
    logsums = -cumsum(log(1-pvalues))
    stepup!(logsums, forwardstop_multiplier, n)
    max(min(logsums, 1.), 0.)
end

forwardstop_multiplier(i::Int, n::Int) = 1/(n-i)
