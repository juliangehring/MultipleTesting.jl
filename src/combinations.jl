### Combination methods for p-values ###

"""
    combine(PValues, <:PValueCombination)

Combine p-values


# Examples

```jldoctest
julia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);

julia> combine(pvals, Fisher())
0.007616871850449092

julia> combine(pvals, Stouffer())
0.00709832618126593
```


# See also

`PValueCombination`s:

[`Fisher`](@ref)
[`Logit`](@ref)
[`Stouffer`](@ref)
[`Tippett`](@ref)
[`Simes`](@ref)
[`Wilkinson`](@ref)
[`Minimum`](@ref)

"""
function combine end

function combine(pValues::AbstractVector{T}, method::M)::T where {T <: AbstractFloat,M <: PValueCombination}
    return combine(PValues(pValues), method)
end


## Fisher combination ##

"""
Fisher's p-value combination


# Examples

```jldoctest
julia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);

julia> combine(pvals, Fisher())
0.007616871850449092
```


# References

Fisher, R.A. (1925).
Statistical methods for research workers
(Genesis Publishing Pvt Ltd).

"""
struct Fisher <: PValueCombination
end

function combine(pValues::PValues{T}, method::Fisher)::T where T <: AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    if minimum(pValues) == 0
        return NaN
    end
    x = -2 * sum(log.(pValues))
    p = ccdf(Chisq(2n), Float64(x))
    return p
end


## Logit combination ##

"""
Logit p-value combination


# Examples

```jldoctest
julia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);

julia> combine(pvals, Logit())
0.006434494635148462
```


# References

Mudholkar, G.S., and George, E.O. (1977).
The Logit Statistic for Combining Probabilities - An Overview
(Rochester University NY, Dept of Statistics).

"""
struct Logit <: PValueCombination
end

function combine(pValues::PValues{T}, method::Logit)::T where T <: AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    if minimum(pValues) == 0 || maximum(pValues) == 1
        return NaN
    end
    c = sqrt((5n + 2) * n * pi^2 / ((5n + 4) * 3))
    x = -sum(log.(pValues ./ (1 .- pValues))) / c
    p = ccdf(TDist(5n + 4), x)
    return p
end


## Stouffer combination ##

"""
Stouffer's p-value combination


# Examples

```jldoctest
julia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);

julia> combine(pvals, Stouffer())
0.00709832618126593

julia> weights = [1.0, 2.0, 0.4, 1.5];

julia> combine(pvals, weights, Stouffer())
0.007331653763696742
```


# References

Stouffer, S.A. (1949).
The American soldier. Vol. 1: Adjustment during army life
(Princeton University Press).

Liptak, T. (1958).
On the combination of independent tests.
Magyar Tud Akad Mat Kutato Int Kozl 3, 171–197.

"""
struct Stouffer <: PValueCombination
end

function combine(pValues::PValues{T}, method::Stouffer)::T where T <: AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    if minimum(pValues) == 0 || maximum(pValues) == 1
        return NaN
    end
    z = cquantile.(Ref(Normal()), pValues)
    z = sum(z) ./ sqrt(n)
    p = ccdf(Normal(), z)
    return p
end

function combine(pValues::PValues{T}, weights::AbstractVector{R}, method::Stouffer)::T where {T <: AbstractFloat,R <: Real}
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    if minimum(pValues) == 0 || maximum(pValues) == 1
        return NaN
    end
    z = cquantile.(Ref(Normal()), pValues) .* weights
    z = sum(z) ./ sqrt(sum(abs2, weights))
    p = ccdf(Normal(), z)
    return p
end

function combine(pValues::AbstractVector{T}, weights::AbstractVector{R}, method::Stouffer)::T where {T <: AbstractFloat,R <: Real}
    return combine(PValues(pValues), weights, method)
end


## Tippett combination ##

"""
Tippett's p-value combination


# Examples

```jldoctest
julia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);

julia> combine(pvals, Tippett())
0.039403990000000055
```


# References

Tippett, L.H.C. (1931). The Methods of Statistics. An introduction mainly for
workers in the biological sciences.

"""
struct Tippett <: PValueCombination
end

function combine(pValues::PValues{T}, method::Tippett)::T where T <: AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    p = 1 - (1 - minimum(pValues))^n
    return p
end


## Simes combination ##

"""
Simes's p-value combination


# Examples

```jldoctest
julia> pvals = PValues([0.01, 0.02, 0.3, 0.5]);

julia> combine(pvals, Simes())
0.04
```


# References

Simes, R.J. (1986). An improved Bonferroni procedure for multiple tests of
significance. Biometrika 73, 751–754.

"""
struct Simes <: PValueCombination
end

function combine(pValues::PValues{T}, method::Simes)::T where T <: AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pValues = sort(pValues)
    p = n * minimum(pValues ./ (1:n))
    return p
end


## Wilkinson combination ##

"""
Wilkinson's p-value combination


# Examples

```jldoctest
julia> pv = PValues([0.01, 0.02, 0.3, 0.5]);

julia> combine(pv, Wilkinson(1))  # combination with rank 1
0.03940399000000003

julia> combine(pv, Wilkinson(4))  # combination with rank 4
0.0625
```


# References

Wilkinson, B. (1951). A statistical consideration in psychological research.
Psychological Bulletin 48, 156.

"""
struct Wilkinson <: PValueCombination
    rank::Integer

    function Wilkinson(rank)
        if rank < 1
            throw(ArgumentError("Rank must be positive."))
        end
        return new(rank)
    end
end

function combine(pValues::PValues{T}, method::Wilkinson)::T where T <: AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    rank = method.rank
    if rank < 1 || rank > n
        throw(ArgumentError("Rank must be in 1,..,$(n)"))
    end
    p_rank = partialsort(pValues, rank)
    p = cdf(Beta(rank, n - rank + 1), Float64(p_rank))
    return p
end


## Generalised minimum combination ##

"""
Minimum of adjusted p-value combination


# Examples

```jldoctest
julia> pv = PValues([0.01, 0.02, 0.3, 0.5]);

julia> combine(pv, Minimum(BenjaminiHochberg()))
0.04

julia> combine(pv, Minimum(ForwardStop()))
0.01005033585350145
```

"""
struct Minimum <: PValueCombination
    adjustment::PValueAdjustment
end

function combine(pValues::PValues{T}, method::Minimum)::T where T <: AbstractFloat
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    padj = adjust(pValues, method.adjustment)
    p = minimum(padj)
    return p
end
