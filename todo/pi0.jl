function smooth_pi0{T<:FloatingPoint}(pValues::Vector{T}, lambda = 0.05:0.05:0.95)
    validPvalues(pValues)
    #validPvalues(lambda) ## TODO check bounds
    n = length(pValues)
    ## CHCK: Only for smoothing (due to cubic spline)?
    if length(lambda) < 4
        throw(ArgumentError())
    end
    if !issorted(lambda)
        throw(DomainError())
    end
    pi0 = Float64[mean(pValues .>= l) / (1-l) for l in lambda]
    ##CoordInterpGrid(lambda, pi0, BCnil, InterpCubic) ## broken
end
