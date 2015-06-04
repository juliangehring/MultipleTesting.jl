## Test util methods ##
module Test_utils

using MultipleTesting
using Base.Test

## validPValues ##
println(" ** ", "validPValues")
@test_throws DomainError MultipleTesting.validPValues([-1.])
@test_throws DomainError MultipleTesting.validPValues([2.])

## reorder ##
println(" ** ", "reorder")
x = [1, 5, 4, 2, 4, 3]
no, oo = MultipleTesting.reorder(x)
@test x[no] == sort(x)
@test x[no][oo] == x

end
