__precompile__()

"""
*MultipleTesting* package

* Adjusting p-values for multiple testing
* Estimating the fraction π₀ of tests under the null hypothesis
* Combination of p-values
* Higher criticism
"""
module MultipleTesting

using StatsBase
import StatsBase: fit

using Distributions
import Distributions: estimate

import SpecialFunctions: digamma

import Random: shuffle!

export
    PValues,
    adjust,
    PValueAdjustment,
    Bonferroni,
    BenjaminiHochberg,
    BenjaminiHochbergAdaptive,
    BenjaminiYekutieli,
    BenjaminiLiu,
    Hochberg,
    Holm,
    Hommel,
    Sidak,
    ForwardStop,
    BarberCandes,
    estimate_pi0,
    Pi0Estimator,
    Storey,
    StoreyBootstrap,
    LeastSlope,
    Oracle,
    TwoStep,
    RightBoundary,
    CensoredBUM,
    CensoredBUMFit,
    BUM,
    BUMFit,
    FlatGrenander,
    ConvexDecreasing,
    ConvexDecreasingFit,
    isin,
    fit,
    BetaUniformMixtureModel,
    PValueCombination,
    combine,
    FisherCombination,
    LogitCombination,
    StoufferCombination,
    TippettCombination,
    SimesCombination,
    WilkinsonCombination,
    MinimumCombination,
    estimate,
    HigherCriticismScores,
    HigherCriticismThreshold


include("types.jl")
include("utils.jl")
include("pval-adjustment.jl")
include("pi0-estimators.jl")
include("model.jl")
include("combinations.jl")
include("higher-criticism.jl")

end
