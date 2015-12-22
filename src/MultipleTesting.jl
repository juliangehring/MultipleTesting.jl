module MultipleTesting

using Compat
using StatsBase
import StatsBase.fit

export
    qValues,
    storey_pi0,
    bootstrap_pi0,
    lsl_pi0,
    bonferroni,
    benjamini_hochberg,
    benjamini_yekutieli,
    holm,
    hommel,
    hochberg,
    qValues,
    isotonicregression

include("utils.jl")
include("pi0-estimators.jl")
include("qvalue.jl")
include("pval-adjustment.jl")

end # module
