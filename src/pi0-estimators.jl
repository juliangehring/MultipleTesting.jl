### estimators for π₀ (pi0)

function estimate_pi0(pValues::AbstractVector{T}, method::M) where {T<:AbstractFloat, M<:Pi0Estimator}
    estimate_pi0(PValues(pValues), method)
end


## Storey estimator

"""
Storey π₀ estimator

**Parameters**

- λ : tuning parameter, FloatingPoint, default: 0.1

**Examples**

```julia
Storey()
Storey(0.1)
```

**References**

Storey, JD (2002). "A Direct Approach to False Discovery Rates." Journal of the
Royal Statistical Society, doi:10.1111/1467-9868.00346

"""
struct Storey <: Pi0Estimator
    λ::Float64

    Storey(λ) = isin(λ, 0, 1) ? new(λ) : throw(DomainError())
end

Storey() = Storey(0.1)

function estimate_pi0(pValues::PValues{T}, pi0estimator::Storey) where T<:AbstractFloat
    lambda = pi0estimator.λ
    pi0 = (sum(pValues .>= lambda) / length(pValues)) / (1 - lambda)
    pi0 = min.(pi0, 1)
    return pi0
end


## Storey bootstrap estimator
"""
Storey closed-form bootstrap π₀ estimator

StoreyBootstrap(λseq, q)

Reference: David Robinson, 2012
"""
struct StoreyBootstrap <: Pi0Estimator
    λseq::Vector{Float64}
    q   ::Float64

    StoreyBootstrap(λseq, q) =
        isin(λseq, 0, 1) && isin(q, 0, 1) ? new(λseq, q) : throw(DomainError())
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
    pi0 = min.(pi0[indmin(mse)], 1)
    pi0
end


## Least SLope (LSL) estimator

"""
Least SLope (LSL) π₀ estimator

LeastSlope()
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
    pi0 = min.( 1/sx + 1, n ) / n
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
    idx = findfirst(d) + 1
    pi0 = min.( 1/s[idx] + 1, n ) / n
    return pi0
end


## Oracle
"""
Oracle π₀

Oracle(π₀)
"""
struct Oracle <: Pi0Estimator
    π0::Float64

    Oracle(π0) = isin(π0, 0, 1) ? new(π0) : throw(DomainError())
end

Oracle() = Oracle(1.0)

function estimate_pi0(pValues::PValues{T}, pi0estimator::Oracle) where T<:AbstractFloat
    pi0estimator.π0
end


## Two-Step estimator: Benjamini, Krieger and Yekutieli (2006)

"""
Two step π₀ estimator

TwoStep(α)

Reference: Benjamini, Krieger and Yekutieli, 2006
"""
struct TwoStep <: Pi0Estimator
    α::Float64
    method::PValueAdjustment

    TwoStep(α, method) = isin(α, 0, 1) ? new(α, method) : throw(DomainError())
end

TwoStep() = TwoStep(0.05)

TwoStep(α) = TwoStep(α, BenjaminiHochberg())

function estimate_pi0(pValues::PValues{T}, pi0estimator::TwoStep) where T<:AbstractFloat
    alpha = pi0estimator.α
    padj = adjust(pValues, pi0estimator.method)
    pi0 = sum(padj .>= (alpha/(1+alpha))) / length(padj)
    return pi0
end


# RightBoundary procedure as defined in Definition 2 of Liang and Nettleton 2012
# "Adaptive and dynamic adaptive procedures for false discovery rate control and estimation"
"""
Right boundary π₀ estimator

RightBoundary(λseq)
"""
struct RightBoundary <: Pi0Estimator
    λseq::Vector{Float64}

    RightBoundary(λseq) =
        isin(λseq, 0, 1) ? new(λseq) : throw(DomainError())
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
    pi0 = pi0_estimates[findfirst(pi0_decrease, true) + 1]
    return min.(pi0, 1)
end


## Censored BUM
"""
Censored BUM π₀ estimator

CensoredBUM(γ0, λ, xtol, maxiter)
"""
struct CensoredBUM <: Pi0Estimator
    γ0::Float64
    λ::Float64
    xtol::Float64
    maxiter::Int64

    function CensoredBUM(γ0, λ, xtol, maxiter)
        if isin(γ0, 0, 1) && isin(λ, 0, 1) && isin(xtol, 0, 1) && maxiter > 0
            new(γ0, λ, xtol, maxiter)
        else
            throw(DomainError())
        end
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
        # explicitly handle denominator of 0 in julia 0.5: min(x, NaN) == x
        if isnan(α)
           break
        end
        γ = max.(min.(γ, 1), 0)
        α = max.(min.(α, 1), 0)
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
        γ = sum(1 .- z) / n  # TODO simplify
        α = -sum(z[idx_right])
        α = α / ( ll * sum(z[idx_left]) + sum(z[idx_right] .* lpr) )
        xl = (1-γ) * (λ^α)
        z[idx_left] = xl ./ (γ*λ + xl)
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
BUM π₀ estimator

BUM(γ0, xtol, maxiter)
"""
struct BUM <: Pi0Estimator
    γ0::Float64
    xtol::Float64
    maxiter::Int64

    function BUM(γ0, xtol, maxiter)
        if isin(γ0, 0, 1) && isin(xtol, 0, 1)
            new(γ0, xtol, maxiter)
        else
            throw(DomainError())
        end
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


## Longest constant interval in the Grenander estimator: Langaas et al., 2005

"""
Flat Grenander π₀ estimator

FlatGrenander()

Estimates π₀ by the longest constant interval in the Grenander estimator

Reference: Langaas et al., 2005: section 4.3
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

ConvexDecreasing(gridsize, xtol, maxiter)
"""
struct ConvexDecreasing <: Pi0Estimator
    gridsize::Int64
    xtol::Float64
    maxiter::Int64

    function ConvexDecreasing(gridsize, xtol, maxiter)
        if gridsize > 0 && isin(xtol, 0, 1) && maxiter > 0
            new(gridsize, xtol, maxiter)
        else
            throw(DomainError())
        end
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
    return indmax( [theta.^-2 * sum(theta .- p[p .< theta]) for theta in t] )
end

function find_theta(t::Vector{Float64}, p::Vector{Float64}, f_p::Vector{Float64})
    return indmax( [theta.^-2 * sum( (theta .- p) .* (p .< theta) ./ f_p ) for theta in t] )
end

function decide(f_p::Vector{Float64}, f_theta_p::Vector{Float64}, ε::Float64)
    idx = f_p .> 0
    return sum( @. ( f_p[idx] - f_theta_p[idx] ) / ( (1-ε) * f_p[idx] + ε * f_theta_p[idx] ) )
end

function triangular_weighting(x::Vector{T}, mid::T) where T<:AbstractFloat
    return @. 2 / mid^2 * (mid-x) * (x<mid)
end
