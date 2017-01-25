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

export
    adjust,
    PValueAdjustmentMethod,
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
    qValues,
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
    isin,
    fit,
    BetaUniformMixtureModel,
    PValues

include("types.jl")
include("utils.jl")
include("pval-adjustment.jl")
include("pi0-estimators.jl")
include("qvalue.jl")
include("model.jl")

end
