using Documenter
using MultipleTesting

makedocs(
    modules = [MultipleTesting],
    format = :html,
    sitename = "MultipleTesting.jl",
    linkcheck = false
)
