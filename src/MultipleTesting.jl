__precompile__()

module MultipleTesting

using StatsBase

export
    qValues,
    storey_pi0,
    bootstrap_pi0,
    lsl_pi0,
    twostep_pi0,
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
    isin

include("types.jl")
include("utils.jl")
include("pval-adjustment.jl")
include("pi0-estimators.jl")
include("qvalue.jl")

end
