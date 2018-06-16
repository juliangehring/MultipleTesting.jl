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

adjust(pValues::PValues{T}, method::Bonferroni) where T<:AbstractFloat = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues{T}, n::Integer, method::Bonferroni) where T<:AbstractFloat
    k = length(pValues)
    check_number_tests(k, n)
    pAdjusted = clamp.(pValues * n, 0, 1)
    return pAdjusted
end


# Benjamini-Hochberg

struct BenjaminiHochberg <: PValueAdjustment
end

adjust(pValues::PValues{T}, method::BenjaminiHochberg) where T<:AbstractFloat = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues{T}, n::Integer, method::BenjaminiHochberg) where T<:AbstractFloat
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    pAdjusted .*= n ./ (1:k)
    stepup!(pAdjusted)
    pAdjusted = clamp.(pAdjusted[originalOrder], 0, 1)
    return pAdjusted
end


# Benjamini-Hochberg Adaptive

struct BenjaminiHochbergAdaptive <: PValueAdjustment
    pi0estimator::Pi0Estimator
end

## default to BenjaminiHochberg
BenjaminiHochbergAdaptive(π0::AbstractFloat) = BenjaminiHochbergAdaptive(Oracle(π0))

BenjaminiHochbergAdaptive() = BenjaminiHochbergAdaptive(1.0)

adjust(pValues::PValues{T}, method::BenjaminiHochbergAdaptive) where T<:AbstractFloat = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues{T}, n::Integer, method::BenjaminiHochbergAdaptive) where T<:AbstractFloat
    π0 = estimate_pi0(pValues, method.pi0estimator)
    pAdjusted = adjust(pValues, n, BenjaminiHochberg()) * π0
    return pAdjusted
end


# Benjamini-Yekutieli

struct BenjaminiYekutieli <: PValueAdjustment
end

adjust(pValues::PValues{T}, method::BenjaminiYekutieli) where T<:AbstractFloat = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues{T}, n::Integer, method::BenjaminiYekutieli) where T<:AbstractFloat
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    pAdjusted .*= harmonic_number(n) .* n ./ (1:k)
    stepup!(pAdjusted)
    pAdjusted = clamp.(pAdjusted[originalOrder], 0, 1)
    return pAdjusted
end


# Benjamini-Liu

struct BenjaminiLiu <: PValueAdjustment
end

adjust(pValues::PValues{T}, method::BenjaminiLiu) where T<:AbstractFloat = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues{T}, n::Integer, method::BenjaminiLiu) where T<:AbstractFloat
    k = length(pValues)
    check_number_tests(k, n)
    if n <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    # a bit more involved because cutoffs at significance α have the form:
    # P_(i) <= 1- [1 - min(1, m/(m-i+1)α)]^{1/(m-i+1)}
    s = n .- (1:k) .+ 1
    pAdjusted = (1 .- (1 .- pAdjusted) .^ s) .* s ./ n
    stepdown!(pAdjusted)
    pAdjusted = clamp.(pAdjusted[originalOrder], 0, 1)
    return pAdjusted
end


# Hochberg

struct Hochberg <: PValueAdjustment
end

adjust(pValues::PValues{T}, method::Hochberg) where T<:AbstractFloat = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues{T}, n::Integer, method::Hochberg) where T<:AbstractFloat
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    pAdjusted .*= (n .- (1:k) .+ 1)
    stepup!(pAdjusted)
    pAdjusted = clamp.(pAdjusted[originalOrder], 0, 1)
    return pAdjusted
end


# Holm

struct Holm <: PValueAdjustment
end

adjust(pValues::PValues{T}, method::Holm) where T<:AbstractFloat = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues{T}, n::Integer, method::Holm) where T<:AbstractFloat
    k = length(pValues)
    check_number_tests(k, n)
    if n <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    pAdjusted .*= (n .- (1:k) .+ 1)
    stepdown!(pAdjusted)
    pAdjusted = clamp.(pAdjusted[originalOrder], 0, 1)
    return pAdjusted
end


# Hommel

struct Hommel <: PValueAdjustment
end

adjust(pValues::PValues{T}, method::Hommel) where T<:AbstractFloat = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues{T}, n::Integer, method::Hommel) where T<:AbstractFloat
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = vcat(pValues[sortedOrder], fill(one(T), n-k))
    lower = n * minimum(pAdjusted./(1:n))
    q = fill(lower, n)
    pa = fill(lower, n)
    for j in (n-1):-1:2
        idx_left = 1:(n-j+1)
        idx_right = (n-j+2):n
        q_right = minimum(view(pAdjusted, idx_right)./(2:j))
        q[idx_left] .= j .* min.(view(pAdjusted, idx_left), q_right)
        q[idx_right] .= q[n-j+1]
        pa .= max.(pa, q)
    end
    pAdjusted = max.(pa[originalOrder], pValues)
    return pAdjusted
end


# Sidak

struct Sidak <: PValueAdjustment
end

adjust(pValues::PValues{T}, method::Sidak) where T<:AbstractFloat = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues{T}, n::Integer, method::Sidak) where T<:AbstractFloat
    check_number_tests(length(pValues), n)
    pAdjusted = clamp.(1 .- (1 .- pValues).^n, 0, 1)
    return pAdjusted
end


# Forward Stop

struct ForwardStop <: PValueAdjustment
end

adjust(pValues::PValues{T}, method::ForwardStop) where T<:AbstractFloat = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues{T}, n::Integer, method::ForwardStop) where T<:AbstractFloat
    k = length(pValues)
    check_number_tests(k, n)
    sortedOrder, originalOrder = reorder(pValues)
    logsums = -cumsum(log.(1 .- pValues[sortedOrder]))
    logsums ./= (1:k)
    stepup!(logsums)
    pAdjusted = clamp.(logsums[originalOrder], 0, 1)
    return pAdjusted
end


# Barber-Candès
struct BarberCandes <: PValueAdjustment
end

function adjust(pValues::PValues{T}, method::BarberCandes) where T<:AbstractFloat
    n = length(pValues)
    # special cases unlike other p-adjust methods
    if n <= 1
        return fill(1.0, size(pValues))
    end
    if maximum(pValues) < 0.5
        return fill(1/n, size(pValues))
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

    stepup!(estimated_fdrs)
    pAdjusted = clamp.(estimated_fdrs[original_order], 0, 1)
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
    stepup!(estimated_fdrs)
    pAdjusted = clamp.(estimated_fdrs[original_order], 0, 1)
    return pAdjusted
end


## internal ##

# step-up / step-down

function stepup!(values::AbstractVector{T}) where T<:AbstractFloat
    accumulate!(min, values, reverse!(values))
    reverse!(values)
    return values
end

function stepdown!(values::AbstractVector{T}) where T<:AbstractFloat
    accumulate!(max, values, values)
    return values
end


function check_number_tests(k::Integer, n::Integer)
    if k > n
        msg = "Number of total tests ($n) is smaller than number of p-values ($k)"
        throw(ArgumentError(msg))
    end
end


harmonic_number(n::Integer) =  digamma(n+1) + γ
