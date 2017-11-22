## abstract types ##

abstract type Pi0Estimator end

abstract type Pi0Fit end

abstract type PValueAdjustment end

abstract type PValueCombination end


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
