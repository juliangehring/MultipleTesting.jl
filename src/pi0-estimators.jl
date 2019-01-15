### estimators for π₀

"""
    estimate_pi0(PValues, Pi0Estimator)

Estimate π₀, the fraction of tests under the null hypothesis


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, StoreyBootstrap())
0.0
julia> estimate_pi0(pvals, FlatGrenander())
0.42553191489361697
```


# See also

`Pi0Estimator`s:

[`Storey`](@ref)
[`StoreyBootstrap`](@ref)
[`LeastSlope`](@ref)
[`Oracle`](@ref)
[`TwoStep`](@ref)
[`RightBoundary`](@ref)
[`CensoredBUM`](@ref)
[`BUM`](@ref)
[`FlatGrenander`](@ref)
[`ConvexDecreasing`](@ref)

"""
function estimate_pi0 end

function estimate_pi0(pValues::AbstractVector{T}, method::M) where {T<:AbstractFloat, M<:Pi0Estimator}
    estimate_pi0(PValues(pValues), method)
end


## Storey estimator

"""
Storey's π₀ estimator


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, Storey())
0.22222222222222224

julia> estimate_pi0(pvals, Storey(0.4))
0.33333333333333337
```


# References

Storey, J.D., Taylor, J.E., and Siegmund, D. (2004). Strong control,
conservative point estimation and simultaneous conservative consistency of false
discovery rates: a unified approach. Journal of the Royal Statistical Society:
Series B (Statistical Methodology) 66, 187–205.

"""
struct Storey <: Pi0Estimator
    λ::Float64

    function Storey(λ)
        isin(λ, 0, 1) || throw(DomainError("λ must be in [0, 1]"))
        return new(λ)
    end
end

Storey() = Storey(0.1)

function estimate_pi0(pValues::PValues{T}, pi0estimator::Storey) where T<:AbstractFloat
    lambda = pi0estimator.λ
    pi0 = (sum(pValues .>= lambda) / length(pValues)) / (1 - lambda)
    pi0 = clamp(pi0, 0, 1)
    return pi0
end


## Storey bootstrap estimator

"""
Storey's closed-form bootstrap π₀ estimator


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, StoreyBootstrap())
0.0

julia> estimate_pi0(pvals, StoreyBootstrap(0.1:0.1:0.9, 0.2))
0.0
```


# References

Robinson, D. (2016). Original Procedure for Choosing λ.
http://varianceexplained.org/files/pi0boot.pdf

Storey, J.D., Taylor, J.E., and Siegmund, D. (2004). Strong control,
conservative point estimation and simultaneous conservative consistency of false
discovery rates: a unified approach. Journal of the Royal Statistical Society:
Series B (Statistical Methodology) 66, 187–205.

"""
struct StoreyBootstrap <: Pi0Estimator
    λseq::Vector{Float64}
    q   ::Float64

    function StoreyBootstrap(λseq, q)
        isin(λseq, 0, 1) || throw(DomainError("λseq must be in [0, 1]"))
        isin(q, 0, 1) || throw(DomainError("q must be in [0, 1]"))
        return new(λseq, q)
    end
end

StoreyBootstrap() = StoreyBootstrap(0.05:0.05:0.95, 0.1)

function estimate_pi0(pValues::PValues{T}, pi0estimator::StoreyBootstrap) where T<:AbstractFloat
    lambdas = pi0estimator.λseq
    q = pi0estimator.q
    n = length(pValues)
    w = [sum(pValues .>= l) for l in lambdas]  # TODO: check if >= or >
    pi0 = w ./ n ./ (1 .- lambdas)
    min_pi0 = quantile(pi0, q)
    mse = (w ./ (n.^2 .* (1 .- lambdas).^2 )) .* (1 .- w/n) + (pi0 .- min_pi0).^2
    pi0 = clamp(pi0[argmin(mse)], 0, 1)
    return pi0
end


## Least SLope (LSL) estimator

"""
Least Slope (LSL) π₀ estimator


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, LeastSlope())
1.0
```


# References

Benjamini, Y., and Hochberg, Y. (2000). On the Adaptive Control of the False
Discovery Rate in Multiple Testing With Independent Statistics. Journal of
Educational and Behavioral Statistics 25, 60–83.

"""
struct LeastSlope <: Pi0Estimator
end

function estimate_pi0(pValues::PValues{T}, pi0estimator::LeastSlope) where T<:AbstractFloat
    n = length(pValues)
    pValues = sort_if_needed(pValues)
    s0 = lsl_slope(1, n, pValues)
    sx = 0
    for i in 2:n
        s1 = lsl_slope(i, n, pValues)
        if (s1 - s0) < 0
            sx = s1
            break
        end
        s0 = s1
    end
    pi0 = min( 1/sx + 1, n ) / n
    return pi0
end

function lsl_slope(i::Integer, n::Integer, pValue::AbstractVector{T}) where T<:AbstractFloat
    s = (1 - pValue[i]) / (n - i + 1)
    return s
end

# alternative, vectorized version
# used for comparison and compactness
function lsl_pi0_vec(pValues::AbstractVector{T}) where T<:AbstractFloat
    n = length(pValues)
    pValues = sort_if_needed(pValues)
    s = (1 .- pValues) ./ (n:-1:1)
    d = diff(s) .< 0
    idx = something(findfirst(d), 0) + 1
    pi0 = min( 1/s[idx] + 1, n ) / n
    return pi0
end


## Oracle

"""
Oracle π₀


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, Oracle(0.5)) # a bit boring...
0.5
```
"""
struct Oracle <: Pi0Estimator
    π0::Float64

    function Oracle(π0)
        isin(π0, 0, 1) || throw(DomainError("π0 must be in [0, 1]"))
        return new(π0)
    end
end

Oracle() = Oracle(1.0)

function estimate_pi0(pValues::PValues{T}, pi0estimator::Oracle) where T<:AbstractFloat
    pi0estimator.π0
end


## Two-Step estimator: Benjamini, Krieger and Yekutieli (2006)

"""
Two-step π₀ estimator


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, TwoStep())
0.2

julia> estimate_pi0(pvals, TwoStep(0.05, BenjaminiLiu()))
0.2

```


# References

Benjamini, Y., Krieger, A.M., and Yekutieli, D. (2006). Adaptive linear step-up
procedures that control the false discovery rate. Biometrika 93, 491–507.

"""
struct TwoStep <: Pi0Estimator
    α::Float64
    adjustment::PValueAdjustment

    function TwoStep(α, method)
        isin(α, 0, 1) || throw(DomainError("α must be in [0, 1]"))
        return new(α, method)
    end
end

TwoStep() = TwoStep(0.05)

TwoStep(α) = TwoStep(α, BenjaminiHochberg())

function estimate_pi0(pValues::PValues{T}, pi0estimator::TwoStep) where T<:AbstractFloat
    alpha = pi0estimator.α
    padj = adjust(pValues, pi0estimator.adjustment)
    pi0 = sum(padj .>= (alpha/(1+alpha))) / length(padj)
    return pi0
end


# RightBoundary procedure as defined in Definition 2 of Liang and Nettleton, 2012

"""
Right boundary π₀ estimator


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, RightBoundary())
0.2127659574468085

julia> estimate_pi0(pvals, RightBoundary(0.1:0.1:0.9))
0.25
```


# References

Liang, K., and Nettleton, D. (2012). Adaptive and dynamic adaptive procedures
for false discovery rate control and estimation. Journal of the Royal
Statistical Society: Series B (Statistical Methodology) 74, 163–182.

"""
struct RightBoundary <: Pi0Estimator
    λseq::Vector{Float64}

    function RightBoundary(λseq)
        isin(λseq, 0, 1) || throw(DomainError("λseq must be in [0, 1]"))
        return new(λseq)
    end
end

# λseq used in Liang, Nettleton 2012
RightBoundary() = RightBoundary([0.02:0.02:0.1; 0.15:0.05:0.95])

function estimate_pi0(pValues::PValues{T}, pi0estimator::RightBoundary) where T<:AbstractFloat
    n = length(pValues)
    λseq = sort_if_needed(pi0estimator.λseq)
    # make sure we catch p-values equal to 1 despite left closure
    # use closed=:left because we have been using >= convention in this package
    # note that original paper uses > convention.
    h = fit(Histogram, pValues, [λseq; Inf], closed=:left)
    pi0_estimates = reverse(cumsum(reverse(h.weights)))./(1 .- λseq)./n
    pi0_decrease = diff(pi0_estimates) .>= 0
    pi0_decrease[end] = true
    pi0 = pi0_estimates[something(findfirst(pi0_decrease), 0) + 1]
    pi0 = clamp(pi0, 0, 1)
    return pi0
end


## Censored BUM

"""
Censored Beta-Uniform Mixture (censored BUM) π₀ estimator


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, CensoredBUM())
0.21052495526400936

```


# References

Markitsis, A., and Lai, Y. (2010). A censored beta mixture model for the
estimation of the proportion of non-differentially expressed genes.
Bioinformatics 26, 640–646.

"""
struct CensoredBUM <: Pi0Estimator
    γ0::Float64
    λ::Float64
    xtol::Float64
    maxiter::Int64

    function CensoredBUM(γ0, λ, xtol, maxiter)
        isin(γ0, 0, 1) || throw(DomainError("γ0 must be in [0, 1]"))
        isin(λ, 0, 1) || throw(DomainError("λ must be in [0, 1]"))
        isin(xtol, 0, 1) || throw(DomainError("xtol must be in [0, 1]"))
        maxiter > 0 || throw(DomainError("maxiter must be a positive number"))
        return new(γ0, λ, xtol, maxiter)
    end
end

CensoredBUM() = CensoredBUM(0.5, 0.05, 1e-6, 10000)
CensoredBUM(γ0, λ) = CensoredBUM(γ0, λ, 1e-6, 10000)

struct CensoredBUMFit <: Pi0Fit
    π0::Float64
    param::Vector{Float64}
    is_converged::Bool
end

function fit(pi0estimator::CensoredBUM, pValues::AbstractVector{T}; kw...) where T<:AbstractFloat
    π0, param, is_converged = cbum_pi0(pValues, pi0estimator.γ0, pi0estimator.λ,
                                       pi0estimator.xtol, pi0estimator.maxiter;
                                       kw...)
    return CensoredBUMFit(π0, param, is_converged)
end

function estimate_pi0(pValues::PValues{T}, pi0estimator::CensoredBUM) where T<:AbstractFloat
    estimate_pi0(fit(pi0estimator, pValues))
end

function estimate_pi0(pi0fit::CensoredBUMFit)
    pi0 = pi0fit.is_converged ? pi0fit.π0 : NaN
    return pi0
end

function cbum_pi0(pValues::AbstractVector{T},
                  γ0::AbstractFloat = 0.5, λ::AbstractFloat = 0.05,
                  xtol::AbstractFloat = 1e-6, maxiter::Integer = 10000) where T<:AbstractFloat
    n = length(pValues)
    idx_right = pValues .>= λ
    n2 = sum(idx_right)
    n1 = n - n2
    sz = (1-γ0)*n
    szr = (1-γ0)*n2
    szl = sz - szr
    pr = pValues[idx_right]
    lpr = log.(pr)
    ll = log(λ)
    zr = fill(1-γ0, n2)
    pi0_old = γ0 = α = γ = Inf
    for i in 1:maxiter
        γ = 1 - sz/n
        α = -szr / ( ll * szl + sum(zr .* lpr) )
        γ = clamp(γ, 0, 1)
        α = clamp(α, 0, 1)
        xl = (1-γ) * (λ^α)
        szl = (xl ./ (γ*λ + xl)) * n1
        xr = (1-γ) * α * pr.^(α-1)
        zr = xr ./ (γ .+ xr)
        szr = sum(zr)
        sz = szl + szr
        pi0_new = γ + (1-γ)*α
        if abs(pi0_new - pi0_old) <= xtol
            return pi0_new, [γ, α], true
        end
        γ0 = γ
        pi0_old = pi0_new
    end
    return NaN, [γ, α], false
end

function cbum_pi0_naive(pValues::AbstractVector{T},
                        γ0::AbstractFloat = 0.5, λ::AbstractFloat = 0.05,
                        xtol::AbstractFloat = 1e-6, maxiter::Integer = 10000) where T<:AbstractFloat
    n = length(pValues)
    z = fill(1-γ0, n)
    idx_left = pValues .< λ
    idx_right = .!idx_left
    pi0_old = γ0 = α = γ = Inf
    # compute constant values only once
    lpr = log.(pValues[idx_right])
    ll = log(λ)
    for i in 1:maxiter
        γ = sum(1 .- z) / n
        α = -sum(z[idx_right])
        α = α / ( ll * sum(z[idx_left]) + sum(z[idx_right] .* lpr) )
        xl = (1-γ) * (λ^α)
        z[idx_left] .= xl / (γ*λ + xl)
        xr = (1-γ) * α * pValues[idx_right].^(α-1)
        z[idx_right] = xr ./ (γ .+ xr)
        pi0_new = γ + (1-γ)*α
        if abs(pi0_new - pi0_old) <= xtol
            return pi0_new, [γ, α], true
        end
        γ0 = γ
        pi0_old = pi0_new
    end
    return NaN, [γ, α], false
end


## BUM

"""
Beta-Uniform Mixture (BUM) π₀ estimator


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, BUM())
0.22802795505154264
```


# References

Pounds, S., and Morris, S.W. (2003). Estimating the occurrence of false
positives and false negatives in microarray studies by approximating and
partitioning the empirical distribution of p-values. Bioinformatics 19,
1236–1242.

"""
struct BUM <: Pi0Estimator
    γ0::Float64
    xtol::Float64
    maxiter::Int64

    function BUM(γ0, xtol, maxiter)
        isin(γ0, 0, 1) || throw(DomainError("γ0 must be in [0, 1]"))
        isin(xtol, 0, 1) || throw(DomainError("xtol must be in [0, 1]"))
        maxiter > 0 || throw(DomainError("maxiter must be a positive number"))
        return new(γ0, xtol, maxiter)
    end
end

BUM() = BUM(0.5, 1e-6, 10000)

BUM(y0::Float64) = BUM(y0, 1e-6, 10000)

struct BUMFit <: Pi0Fit
    π0::Float64
    param::Vector{Float64}
    is_converged::Bool
end

function fit(pi0estimator::BUM, pValues::AbstractVector{T}; kw...) where T<:AbstractFloat
    π0, param, is_converged = cbum_pi0(pValues, pi0estimator.γ0, eps(),
                                       pi0estimator.xtol, pi0estimator.maxiter;
                                       kw...)
    return BUMFit(π0, param, is_converged)
end

function estimate_pi0(pValues::PValues{T}, pi0estimator::BUM) where T<:AbstractFloat
    estimate_pi0(fit(pi0estimator, pValues))
end

function estimate_pi0(pi0fit::BUMFit)
    pi0 = pi0fit.is_converged ? pi0fit.π0 : NaN
    return pi0
end


## Longest constant interval in the Grenander estimator

"""
Flat Grenander π₀ estimator

Estimates π₀ by finding the longest constant interval in the Grenander estimator.


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, FlatGrenander())
0.42553191489361697
```


# References

Langaas, M., Lindqvist, B.H., and Ferkingstad, E. (2005). Estimating the
proportion of true null hypotheses, with application to DNA microarray data.
Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67,
555–572.

"""
struct FlatGrenander <: Pi0Estimator
end

function estimate_pi0(pValues::PValues{T}, pi0estimator::FlatGrenander) where T<:AbstractFloat
    p, f, F = grenander(pValues)
    pi0 = longest_constant_interval(p, f)
    return pi0
end

function longest_constant_interval(p::AbstractVector{T}, f::AbstractVector{T}) where T<:AbstractFloat
    p = [0.0; p]
    f = [Inf; f]
    i2 = length(f)
    Δp_max = Δp = 0.0
    pi0 = 1
    for i1 in (length(f)-1):-1:1
        if f[i2] ≈ f[i1] # within constant interval
            Δp = p[i2] - p[i1]
        else
            if Δp >= Δp_max
                Δp_max = Δp
                pi0 = f[i2]
            end
            i2 = i1
            Δp = 0.0
        end
        if f[i1] > 1
            break
        end
    end
    return pi0
end


## Convex Decreasing π₀ estimator

"""
Convex Decreasing π₀ estimator


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.002, 0.01, 0.03, 0.5]);

julia> estimate_pi0(pvals, ConvexDecreasing())
0.013007051336745304
```


# References

Langaas, M., Lindqvist, B.H., and Ferkingstad, E. (2005). Estimating the
proportion of true null hypotheses, with application to DNA microarray data.
Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67,
555–572.

"""
struct ConvexDecreasing <: Pi0Estimator
    gridsize::Int64
    xtol::Float64
    maxiter::Int64

    function ConvexDecreasing(gridsize, xtol, maxiter)
        gridsize > 0 || throw(DomainError("gridsize must be a positive number"))
        isin(xtol, 0, 1) || throw(DomainError("xtol must be in [0, 1]"))
        maxiter > 0 || throw(DomainError("maxiter must be a positive number"))
        return new(gridsize, xtol, maxiter)
    end
end

ConvexDecreasing() = ConvexDecreasing(100, 1e-5, 10000)

struct ConvexDecreasingFit <: Pi0Fit
    π0::Float64
    param::Vector{Float64}
    is_converged::Bool
end

function fit(pi0estimator::ConvexDecreasing, pValues::AbstractVector{T}) where T<:AbstractFloat
    π0, param, is_converged = convex_decreasing(pValues,
                                                pi0estimator.gridsize,
                                                pi0estimator.xtol,
                                                pi0estimator.maxiter)
    return ConvexDecreasingFit(π0, param, is_converged)
end

function estimate_pi0(pValues::PValues{T}, pi0estimator::ConvexDecreasing) where T<:AbstractFloat
    estimate_pi0(fit(pi0estimator, pValues))
end

function estimate_pi0(pi0fit::ConvexDecreasingFit)
    pi0 = pi0fit.is_converged ? pi0fit.π0 : NaN
    return pi0
end

function convex_decreasing(pValues::AbstractVector{T},
                           gridsize::Integer = 100,
                           xtol::AbstractFloat = 1e-5,
                           maxiter::Integer = 10000) where T<:AbstractFloat

    n = length(pValues)
    p = sort(pValues)
    dx = 1 / gridsize
    t = collect(dx:dx:1)
    x = collect(0:dx:1)
    f = fill(one(T), gridsize+1)
    f_p = fill(one(T), n)
    theta = dx * find_theta(t, p)
    f_theta = triangular_weighting(x, theta)
    f_theta_p = triangular_weighting(p, theta)
    idx_lower = round.(Int, round.(gridsize .* p, RoundDown) .+ 1)
    p_upper = round.(gridsize .* p, RoundUp)/gridsize
    idx_upper = round.(Int, gridsize .* p_upper) .+ 1
    px = p_upper .- p
    thetas = T[]
    pi0_new = pi0_old = Inf
    for j in 1:maxiter
        if sum((f_p.-f_theta_p)./f_p) > 0
            ε = 0.0
        else
            l = 0.0
            u = 1
            while abs(u-l) > xtol
                ε = (l+u)/2.0
                if decide(f_p, f_theta_p, ε) < 0
                    l = ε
                else
                    u = ε
                end
            end
        end
        @. f = (1-ε)*f + ε*f_theta
        @. f_p = (1-ε)*f_p + ε*f_theta_p
        theta = dx * find_theta(t, p, f_p)
        f_theta .= triangular_weighting(x, theta)
        pi0_new = f[end]
        if abs(pi0_new - pi0_old) <= xtol
            return pi0_new, thetas, true
        end
        pi0_old = pi0_new
        f_theta_p .= triangular_weighting(p, theta)
        if sum(f_theta_p ./ f_p) < sum(1 ./ f_p)
            theta = 0.0
            f_theta .= fill(1, size(f_theta))
            f_theta_p .= fill(1, size(f_theta_p))
        end
        if !(theta in thetas)
            append!(thetas, theta)
            sort!(thetas)
        end
    end
    return NaN, thetas, false
end

function find_theta(t::Vector{Float64}, p::Vector{Float64})
    return argmax( [theta.^-2 * sum(theta .- p[p .< theta]) for theta in t] )
end

function find_theta(t::Vector{Float64}, p::Vector{Float64}, f_p::Vector{Float64})
    return argmax( [theta.^-2 * sum( (theta .- p) .* (p .< theta) ./ f_p ) for theta in t] )
end

function decide(f_p::Vector{Float64}, f_theta_p::Vector{Float64}, ε::Float64)
    idx = f_p .> 0
    return sum( @. ( f_p[idx] - f_theta_p[idx] ) / ( (1-ε) * f_p[idx] + ε * f_theta_p[idx] ) )
end

function triangular_weighting(x::Vector{T}, mid::T) where T<:AbstractFloat
    return @. 2 / mid^2 * (mid-x) * (x<mid)
end
