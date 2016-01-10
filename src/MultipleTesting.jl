__precompile__()

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
    bonferroni,
    benjamini_hochberg,
    benjamini_yekutieli,
    holm,
    hommel,
    hochberg,
    sidak,
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
