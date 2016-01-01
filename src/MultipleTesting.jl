module MultipleTesting

using StatsBase

export
    qValues,
    storey_pi0,
    bootstrap_pi0,
    lsl_pi0,
    adjust,
    PValueAdjustmentMethod,
    Bonferroni,
    BenjaminiHochberg,
    BenjaminiHochbergAdaptive,
    BenjaminiYekutieli,
    Hochberg,
    Holm,
    Hommel,
    bonferroni,
    benjamini_hochberg,
    benjamini_yekutieli,
    holm,
    hommel,
    hochberg,
    qValues,
    isotonicregression,
    estimate_pi0,
    Pi0Estimator,
    Storey,
    StoreyBootstrap,
    LeastSlope,
    Oracle,
    isin


include("utils.jl")
include("pi0-estimators.jl")
include("qvalue.jl")
include("pval-adjustment.jl")

end # module
