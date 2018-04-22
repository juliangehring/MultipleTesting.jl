## p-value adjustment methods ##

# promotion from float vectors to PValues type

function adjust(pValues::Vector{T}, method::M) where {T<:AbstractFloat, M<:PValueAdjustment}
    adjust(PValues(pValues), method)
end

function adjust(pValues::Vector{T}, n::Integer, method::M) where {T<:AbstractFloat, M<:PValueAdjustment}
    adjust(PValues(pValues), n, method)
end

# Bonferroni

struct Bonferroni <: PValueAdjustment
end

adjust(pValues::PValues, method::Bonferroni) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::Bonferroni)
    k = length(pValues)
    check_number_tests(k, n)
    return min.(pValues * n, 1)
end


# Benjamini-Hochberg

struct BenjaminiHochberg <: PValueAdjustment
end

adjust(pValues::PValues, method::BenjaminiHochberg) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::BenjaminiHochberg)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, bejamini_hochberg_step, k, n)
    pAdjusted = min.(sortedPValues[originalOrder], 1)
    return pAdjusted
end

bejamini_hochberg_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p * n/(k-i)


# Benjamini-Hochberg Adaptive

struct BenjaminiHochbergAdaptive <: PValueAdjustment
    pi0estimator::Pi0Estimator
end

## default to BenjaminiHochberg
BenjaminiHochbergAdaptive(π0::AbstractFloat) = BenjaminiHochbergAdaptive(Oracle(π0))

BenjaminiHochbergAdaptive() = BenjaminiHochbergAdaptive(1.0)

adjust(pValues::PValues, method::BenjaminiHochbergAdaptive) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::BenjaminiHochbergAdaptive)
    π0 = estimate_pi0(pValues, method.pi0estimator)
    pAdjusted = adjust(pValues, n, BenjaminiHochberg()) * π0
    return pAdjusted
end


# Benjamini-Yekutieli

struct BenjaminiYekutieli <: PValueAdjustment
end

adjust(pValues::PValues, method::BenjaminiYekutieli) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::BenjaminiYekutieli)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, benjamini_yekutieli_step, k, n)
    pAdjusted = min.(sortedPValues[originalOrder], 1)
    return pAdjusted
end

benjamini_yekutieli_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p * harmonic_number(n) * n/(k-i)


# Benjamini-Liu

struct BenjaminiLiu <: PValueAdjustment
end

adjust(pValues::PValues, method::BenjaminiLiu) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::BenjaminiLiu)
    k = length(pValues)
    check_number_tests(k, n)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepdown!(sortedPValues, benjamini_liu_step, k, n)
    pAdjusted = min.(sortedPValues[originalOrder], 1)
    return pAdjusted
end

function benjamini_liu_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer)
    # a bit more involved because cutoffs at significance α have the form:
    # P_(i) <= 1- [1 - min(1, m/(m-i+1)α)]^{1/(m-i+1)}
    s = n-i+1
    return (1 - (1-p)^s) * s / n
end


# Hochberg

struct Hochberg <: PValueAdjustment
end

adjust(pValues::PValues, method::Hochberg) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::Hochberg)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, hochberg_step, k, n)
    pAdjusted = min.(sortedPValues[originalOrder], 1)
    return pAdjusted
end

hochberg_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p * (n-k+i+1)


# Holm

struct Holm <: PValueAdjustment
end

adjust(pValues::PValues, method::Holm) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::Holm)
    k = length(pValues)
    check_number_tests(k, n)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepdown!(sortedPValues, holm_step, k, n)
    pAdjusted = min.(sortedPValues[originalOrder], 1)
    return pAdjusted
end

holm_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p * (n-i+1)


# Hommel

struct Hommel <: PValueAdjustment
end

adjust(pValues::PValues, method::Hommel) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::Hommel)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    pValues = vcat(pValues, fill(1.0, n-k))  # TODO avoid sorting of ones
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    q = fill(minimum(n .* pValues./(1:n)), n)
    pa = fill(q[1], n)
    for j in (n-1):-1:2
        ij = 1:(n-j+1)
        i2 = (n-j+2):n
        q1 = minimum(j .* sortedPValues[i2]./((2:j)))
        q[ij] = min.(j .* sortedPValues[ij], q1)
        q[i2] = q[n-j+1]
        pa = max.(pa, q)
    end
    pAdjusted = max.(pa, sortedPValues)[originalOrder[1:k]]
    return pAdjusted
end


# Sidak

struct Sidak <: PValueAdjustment
end

adjust(pValues::PValues, method::Sidak) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::Sidak)
    check_number_tests(length(pValues), n)
    pAdjusted = min.(1 .- (1 .- pValues).^n, 1)
    return pAdjusted
end


# Forward Stop

struct ForwardStop <: PValueAdjustment
end

adjust(pValues::PValues, method::ForwardStop) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::ForwardStop)
    k = length(pValues)
    check_number_tests(k, n)
    logsums = -cumsum(log.(1 .- pValues))
    stepup!(logsums, forwardstop_step, k, n)
    pAdjusted = max.(min.(logsums, 1), 0)
    return pAdjusted
end

forwardstop_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p * 1/(k-i)


# Barber-Candès
struct BarberCandes <: PValueAdjustment
end

function adjust(pValues::PValues, method::BarberCandes)
    n = length(pValues)
    if n <= 1
        return fill(1.0, size(pValues)) # unlike other p-adjust methods
    end

    sorted_indexes, original_order = reorder(pValues)
    estimated_fdrs = pValues[sorted_indexes]

    Rt = 1 # current number of discoveries
    Vt = 1 # estimated false discoveries at t

    left_pv = estimated_fdrs[1]
    right_pv = estimated_fdrs[n]

    while left_pv < 0.5
        while (1 - right_pv <= left_pv) && (right_pv >= 0.5) && (Vt + Rt <= n)
            estimated_fdrs[n-Vt+1] = 1.0
            right_pv = estimated_fdrs[n-Vt]
            Vt += 1
        end
        estimated_fdrs[Rt] = Vt/Rt
        Rt += 1
        left_pv = estimated_fdrs[Rt]
    end

    while (right_pv >= 0.5) && (Vt + Rt <= n+1)
      estimated_fdrs[n-Vt+1] = 1.0
      right_pv = (Vt + Rt <= n) ? estimated_fdrs[n-Vt] : 0.0
      Vt += 1
    end

    stepup!(estimated_fdrs, identity_step, n, n) # just monotonize, no multiplier needed
    pAdjusted = min.(estimated_fdrs[original_order], 1)
    return pAdjusted
end

# as test, inefficient implementation
function barber_candes_brute_force(pValues::AbstractVector{T}) where T<:AbstractFloat
    n = length(pValues)
    sorted_indexes, original_order = reorder(pValues)
    sorted_pValues = pValues[sorted_indexes]
    estimated_fdrs = fill(1.0, size(pValues))
    for (i,pv) in enumerate(sorted_pValues)
        if pv >= 0.5
            break
        else
            estimated_fdrs[i] = (sum((1 .- pValues) .<= pv) + 1)/i
        end
    end
    stepup!(estimated_fdrs, identity_step, n, n) # just monotonize, no multiplier needed
    pAdjusted = min.(estimated_fdrs[original_order], 1)
    return pAdjusted
end

identity_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p


## internal ##

# step-up / step-down

function stepup!(sortedPValues::AbstractVector{T}, stepfun::Function, k::Integer, n::Integer) where T<:AbstractFloat
    sortedPValues[k] = stepfun(sortedPValues[k], 0, k, n)
    for i in 1:(k-1)
        sortedPValues[k-i] = min.(sortedPValues[k-i+1], stepfun(sortedPValues[k-i], i, k, n))
    end
    return sortedPValues
end

function stepdown!(sortedPValues::AbstractVector{T}, stepfun::Function, k::Integer, n::Integer) where T<:AbstractFloat
    sortedPValues[1] = stepfun(sortedPValues[1], 1, k, n)
    for i in 2:k
        sortedPValues[i] = max.(sortedPValues[i-1], stepfun(sortedPValues[i], i, k, n))
    end
    return sortedPValues
end


function check_number_tests(k::Integer, n::Integer)
    if k > n
        msg = "Number of total tests ($n) is smaller than number of p-values ($k)"
        throw(ArgumentError(msg))
    end
end


harmonic_number(n::Integer) =  digamma(n+1) + γ
