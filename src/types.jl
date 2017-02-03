## abstract types ##

abstract Pi0Estimator

abstract Pi0Fit

abstract PValueAdjustmentMethod

# consistent naming
abstract PValueCombinationMethod


## concrete types ##

## PValues

type PValues{T<:AbstractFloat} <: AbstractVector{T}
    values::AbstractVector{T}
    min::T
    max::T

    function(::Type{PValues}){T}(values::AbstractVector{T})
        min, max = extrema(values)
        if min < 0.0 || max > 1.0
            throw(DomainError())
        end
        new{T}(values, min, max)
    end
end

Base.convert{T<:AbstractFloat}(::Type{PValues}, x::AbstractVector{T}) = PValues(x)

Base.size(pv::PValues) = (length(pv.values), )
Base.linearindexing{T<:PValues}(::Type{T}) = Base.LinearFast()
Base.getindex(pv::PValues, i::Int) = pv.values[i]

Base.values(pv::PValues) = pv.values
Base.minimum(pv::PValues) = pv.min
Base.maximum(pv::PValues) = pv.max
Base.extrema(pv::PValues) = (minimum(pv), maximum(pv))
Base.length(pv::PValues) = length(pv.values)
