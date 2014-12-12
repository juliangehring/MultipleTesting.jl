module MultipleTesting

using StatsBase

export
    qValues,
    storey_pi0,
    bootstrap_pi0,
    lsl_pi0,
    bonferroni,
    benjamini_hochberg,
    benjamini_yekutieli,
    holm,
    hochberg,
    validPValues

include("utils.jl")
include("pval-adjustment.jl")
include("pi0-estimators.jl")
include("qvalue.jl")

end # module
