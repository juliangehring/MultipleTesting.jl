#__precompile__()

"""
*MultipleTesting* package

* Adjusting p-values for multiple testing
* Estimating the fraction π0 of tests under the null hypothesis
"""
module MultipleTesting

using StatsBase
import StatsBase: fit

export
    storey_pi0,
    bootstrap_pi0,
    lsl_pi0,
    twostep_pi0,
    rightboundary_pi0,
    adjust,
    PValueAdjustmentMethod,
    Bonferroni,
    BenjaminiHochberg,
    BenjaminiHochbergAdaptive,
    BenjaminiYekutieli,
    Hochberg,
    Holm,
    Hommel,
    Sidak,
    ForwardStop,
    bonferroni,
    benjamini_hochberg,
    benjamini_yekutieli,
    holm,
    hommel,
    hochberg,
    sidak,
    forwardstop,
    qValues,
    isotonicregression,
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
    isin,
    fit

include("types.jl")
include("utils.jl")
include("pval-adjustment.jl")
include("pi0-estimators.jl")
include("qvalue.jl")

end
