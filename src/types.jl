## abstract types ##

abstract type Pi0Estimator end

abstract type Pi0Fit end

abstract type PValueAdjustment end

abstract type PValueCombination end

abstract Alternative


## concrete types ##

## PValues

immutable PValues{T<:AbstractFloat} <: AbstractVector{T}
    values::AbstractVector{T}
    min::T
    max::T

    function(::Type{PValues}){T}(values::AbstractVector{T})
        min, max = extrema(values)
        if min < 0.0 || max > 1.0
            throw(DomainError())
        end
        new{T}(copy(values), min, max)
    end
end

Base.convert{T<:AbstractFloat}(::Type{PValues}, x::AbstractVector{T}) = PValues(x)

Base.size(pv::PValues) = (length(pv.values), )
Base.IndexStyle{T<:PValues}(::Type{T}) = IndexLinear()
Base.getindex(pv::PValues, i::Integer) = pv.values[i]
Base.setindex!(pv::PValues, x::AbstractFloat, i::Integer) =
    throw(ErrorException("Modification of values is not permitted"))

Base.values(pv::PValues) = pv.values
Base.minimum(pv::PValues) = pv.min
Base.maximum(pv::PValues) = pv.max
Base.extrema(pv::PValues) = (minimum(pv), maximum(pv))
Base.length(pv::PValues) = length(pv.values)


## ZScores

type ZScores{T<:AbstractFloat} <: AbstractVector{T}
    values::AbstractVector{T}

    function(::Type{ZScores}){T}(values::AbstractVector{T})
        new{T}(values)
    end
end

Base.convert{T<:AbstractFloat}(::Type{ZScores}, x::AbstractVector{T}) = ZScores(x)

Base.size(zs::ZScores) = (length(zs.values), )
Base.linearindexing{T<:ZScores}(::Type{T}) = Base.LinearFast()
Base.getindex(zs::ZScores, i::Int) = zs.values[i]

Base.values(zs::ZScores) = zs.values


# transformation between PValues and ZScores

immutable Upper <: Alternative end

immutable Lower <: Alternative end

immutable Both  <: Alternative end


# PValues to ZScores

function transform(::Type{PValues}, zs::ZScores, alternative::Type{Lower})
    p = cdf(Normal(), zs)
    return PValues(p)
end

function transform(::Type{PValues}, zs::ZScores, alternative::Type{Upper})
    p = ccdf(Normal(), zs)
    return PValues(p)
end

function transform(::Type{PValues}, zs::ZScores, alternative::Type{Both})
    p = 2 * min( ccdf(Normal(), zs), cdf(Normal(), zs) )
    return PValues(p)
end

transform(pv::Type{PValues}, zs::ZScores) = transform(pv, zs, Both)

transform(T::Type{PValues}, zs::ZScores, alt::Alternative) =
    transform(T, zs, typeof(alt))


# ZScores to PValues

function transform(::Type{ZScores}, pv::PValues, alternative::Type{Upper})
    z = cquantile(Normal(), pv)
    return ZScores(z)
end

function transform(::Type{ZScores}, pv::PValues, alternative::Type{Lower})
    z = quantile(Normal(), pv)
    return ZScores(z)
end

transform(T::Type{ZScores}, pv::PValues, alt::Alternative) =
    transform(T, pv, typeof(alt))
