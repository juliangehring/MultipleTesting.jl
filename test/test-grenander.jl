## Test π0 estimator methods ##
module Test_grenander

using MultipleTesting
using Base.Test


println(" ** ", "isotonic regression")
#pooled adjacent violators example page 10 robertson
@test_approx_eq isotonicregression([22.5; 23.333; 20.833; 24.25], [3.0;3.0;3.0;2.0]) [22.222; 22.222; 22.222; 24.25]
#if input already ordered , then output should be the same
@test_approx_eq isotonicregression([1.0;2.0;3.0]) [1.0;2.0;3.0]

end
