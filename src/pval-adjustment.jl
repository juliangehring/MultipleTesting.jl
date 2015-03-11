# add interface for this as well

abstract MultipleTest
abstract FDRTest  <: MultipleTest
abstract FWERTest <: MultipleTest

# define all Multiple Tests we have below

immutable BenjaminiHochbergTest    <: FDRTest  end
immutable BenjaminiYekutieliTest   <: FDRTest  end
immutable HochbergTest             <: FWERTest end
immutable BonferroniTest           <: FWERTest end
immutable HolmTest                 <: FWERTest end
immutable HommelTest               <: FWERTest end

const BenjaminiHochberg  = BenjaminiHochbergTest()
const BenjaminiYekutieli = BenjaminiYekutieliTest()
const Hochberg           = HochbergTest()
const Bonferroni         = BonferroniTest()
const Holm               = HolmTest()
const Hommel             = HommelTest()

# set default to benjamini hochberg
padjust{T<:FloatingPoint}(pvalues::Vector{T}) = padjust(BenjaminiHochberg, pvalues)

padjust{T<:FloatingPoint}(::BenjaminiHochbergTest,  pvalues::Vector{T})= benjamini_hochberg(pvalues)
padjust{T<:FloatingPoint}(::BenjaminiYekutieliTest, pvalues::Vector{T}) = benjamini_yekutieli(pvalues)
padjust{T<:FloatingPoint}(::HochbergTest, pvalues::Vector{T}) = hochberg(pvalues)
padjust{T<:FloatingPoint}(::BonferroniTest, pvalues::Vector{T}) = bonferroni(pvalues)
padjust{T<:FloatingPoint}(::HolmTest, pvalues::Vector{T}) = holm(pvalues)
padjust{T<:FloatingPoint}(::HommelTest, pvalues::Vector{T}) = hommel(pvalues)


# now also introduce weighted methods using multiple dispatch
# a weight vector for multiple testing has the following properties
# sum(ws) = length(ws), ws[i] >=0 \forall i
# also its division should have special properties (0/0=0, x/0 = 1, x/w <= 1)
# therefore we define our own MTPWeightVec type
# code below ugly prototype, should be improved... 
# nicest thing would be if StatsBase would define an AbstractWeightVec, so that
# we can just define MTPWeightVec <: AbstractWeightVec and inherit a lot of the properties


immutable MTPWeight <: Real
  value::Float64
end

immutable MTPWeightVec
  values::Vector{MTPWeight}
end

Base.length(ws::MTPWeightVec) = length(ws.values)
Base.getindex(ws::MTPWeightVec, i) = getindex(ws.values, i)

function MTPWeightVec(ws::Vector{Float64})
  ws_sum = sum(ws)
  m = length(ws)
  # add check if all ws >= 0 and ws_sum > 0
  ws = ws/ws_sum*m
  ws = map((w)->MTPWeight(w), ws)
  MTPWeightVec(ws)
end

function /(pv::Float64, w::MTPWeight)
  if pv == zero(Float64)
    weightedpv = zero(Float64)
  else
    weightedpv = min(one(Float64), pv/w.value)
  end
  weightedpv
end

function /(pv::Vector{Float64}, ws::MTPWeightVec)
  m=length(pv)
  weighted_pv = Array(Float64,m)
  for i=1:m
    weighted_pv[i] = pv[i]/ws[i]
  end
  weighted_pv
end

# make MTP with weights available only for tests for which it has been proven that the resulting procedures still control the FDR/FWER
# For Bonferroni, Holm, BenjaminiHochberg it was shown in:
# "False discovery control with p-value weighting", 2006, CR Genovese, K Roeder, L Wasserman

for weighted_test = (BenjaminiHochbergTest,BonferroniTest,HolmTest)
  @eval begin
    function padjust{T<:FloatingPoint}(tst::($weighted_test),pv::Vector{T},ws::MTPWeightVec)
      padjust(tst, pv/ws)
    end
  end
end
## p-value adjustment methods ##

function bonferroni{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    return min(pValues * length(pValues), 1.)
end

@doc """
# Benjamini-Hochberg p-value adjustement


## Usage

    pAdjusted = benjamini_hochberg{T<:FloatingPoint}(pValues::Vector{T})

Input arguments:

- `pValues`: Vector of p-values that should be adjusted

Return values:

- `pAdjusted`: Vector of adjusted p-values, matching the input `pValues`.


## References

Benjamini, Y. and Hochberg, Y. (1995):
Controlling the false discovery rate: A practical and powerful approach to multiple testing.
Journal of the Royal Statistical Society

http://en.wikipedia.org/wiki/False_discovery_rate#Benjamini.E2.80.93Hochberg_procedure

## Examples

    pOld = rand(100)
    pNew = benjamini_hochberg(pOld)

""" ->
function benjamini_hochberg{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, bejamini_hochberg_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

bejamini_hochberg_multiplier(i::Int, n::Int) = n/(n-i)


function benjamini_hochberg{T<:FloatingPoint}(pValues::Vector{T}, pi0::T)
    validPValues([pi0])
    benjamini_hochberg(pValues) .* pi0
end


function holm{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepdown!(sortedPValues, holm_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

holm_multiplier(i::Int, n::Int) = (n-i+1)


function benjamini_yekutieli{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, benjamini_yekutieli_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

function benjamini_yekutieli_multiplier(i::Int, n::Int)
    c = sum([1/i for i in 1:n])
    return ((n*c)/(n-i))
end


function hochberg{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, hochberg_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

hochberg_multiplier(i::Int, n::Int) = (i+1)


function hommel{T<:FloatingPoint}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    q = fill(minimum(n .* pValues./[1:n; ]), n)
    pa = fill(q[1], n)
    for j in (n-1):-1:2
        ij = 1:(n-j+1)
        i2 = (n-j+2):n
        q1 = minimum(j .* sortedPValues[i2]./([2:j; ]))
        q[ij] = min(j .* sortedPValues[ij], q1)
        q[i2] = q[n-j+1]
        pa = max(pa, q)
    end
    max(pa, sortedPValues)[originalOrder]
end
