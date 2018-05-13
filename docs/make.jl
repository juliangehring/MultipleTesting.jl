using Documenter
using MultipleTesting

makedocs(
    modules = [MultipleTesting],
    format = :html,
    sitename = "MultipleTesting.jl",
    doctest = :fix,
    checkdocs = :exports,
    linkcheck = false,
    pages = [
        "Home" => "index.md",
        "Library" => Any[
            "adjustment.md",
            "combination.md",
            "pi0.md",
            "higher-criticism.md",
            "models.md",
        ]
    ]
)
