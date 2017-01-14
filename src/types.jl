## abstract types ##

abstract Pi0Estimator

abstract Pi0Fit

abstract PValueAdjustmentMethod


## types ##

# PValues

type PValues{T<:AbstractFloat} <: AbstractVector{T}
    values::AbstractVector{T}
    min::T
    max::T

    function PValues(values::AbstractVector{T}, min::T, max::T)
        if min < 0.0 || max > 1.0
            throw(DomainError())
        end
        new(values, min, max)
    end
end

PValues{T<:AbstractFloat}(values::AbstractVector{T}, min::T, max::T) = PValues{T}(values, min, max)

function PValues{T<:AbstractFloat}(values::AbstractVector{T})
    minp, maxp = extrema(values)
    PValues(values, minp, maxp)
end

Base.convert{T<:AbstractFloat}(::Type{PValues}, x::AbstractVector{T}) = PValues(x) # TODO define test

Base.size(pv::PValues) = (length(pv.values), )
Base.linearindexing{T<:PValues}(::Type{T}) = Base.LinearFast()
Base.getindex(pv::PValues, i::Int) = pv.values[i]

Base.values(pv::PValues) = pv.values
Base.minimum(pv::PValues) = pv.min
Base.maximum(pv::PValues) = pv.max
Base.extrema(pv::PValues) = (minimum(pv), maximum(pv))
# not essential: defaults to prod(size(pv))
# store in extra field 'n'?
Base.length(pv::PValues) = length(pv.values)
