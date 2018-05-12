using Documenter
using MultipleTesting

makedocs(
    modules = [MultipleTesting],
    format = :html,
    sitename = "MultipleTesting.jl",
    doctest = :fix
)
