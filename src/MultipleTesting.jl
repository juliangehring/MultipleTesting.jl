__precompile__()

"""
*MultipleTesting* package

* Adjusting p-values for multiple testing
* Estimating the fraction π0 of tests under the null hypothesis
"""
module MultipleTesting

using StatsBase
import StatsBase: fit

using Distributions
import Distributions: estimate

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
    QValues,
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
include("qvalue.jl")
include("model.jl")
include("combinations.jl")
include("higher-criticism.jl")

end
