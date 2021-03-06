## utility functions ##

function valid_pvalues(x::AbstractVector{T}) where T <: AbstractFloat
    if !isin(x)
        throw(DomainError("p-values must all be in [0, 1]"))
    end
end


function isin(x::Real, lower::Real = 0., upper::Real = 1.)
    x >= lower && x <= upper
end

function isin(x::AbstractVector{T}, lower::Real = 0., upper::Real = 1.) where T <: Real
    ex = extrema(x)
    ex[1] >= lower && ex[2] <= upper
end


function isotonic_regression(y::AbstractVector{T}, weights::AbstractVector{T}) where T <: AbstractFloat
    n = length(y)
    if n <= 1
        return y
    end
    if n != length(weights)
        throw(DimensionMismatch("Lengths of values and weights mismatch"))
    end
    @inbounds begin
        n -= 1
        while true
            i = 1
            is_pooled = false
            while i <= n
                k = i
                while k <= n && y[k] >= y[k + 1]
                    k += 1
                end
                if y[i] != y[k]
                    numerator = 0.0
                    denominator = 0.0
                    for j in i:k
                        numerator += y[j] * weights[j]
                        denominator += weights[j]
                    end
                    m = numerator / denominator
                    for j in i:k
                        y[j] = m
                    end
                    is_pooled = true
                end
                i = k + 1
            end
            if !is_pooled
                break
            end
        end
    end
    return y
end

function isotonic_regression(y::AbstractVector{T}) where T <: AbstractFloat
    isotonic_regression(y, fill(1.0, size(y)))
end


function grenander(pv::AbstractVector{T}) where T <: AbstractFloat
    pv_sorted = sort(pv)
    # ecdf that handles duplicated values
    pv_sorted_unique, counts = rle(pv_sorted)
    ecdf_value = cumsum(counts)
    ecdf_value = ecdf_value ./ ecdf_value[end]

    Δx = diff(pv_sorted_unique)
    Δy = diff(ecdf_value)

    f = Δy ./ Δx
    f = -isotonic_regression(-f, Δx)
    F  = ecdf_value[1] .+ vcat(0, cumsum(f .* Δx))
    f = push!(f, f[end])

    return pv_sorted_unique, f, F
end
