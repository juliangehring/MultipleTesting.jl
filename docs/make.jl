using Documenter
using MultipleTesting

makedocs(
    modules = [MultipleTesting],
    format = :html,
    sitename = "MultipleTesting.jl",
    doctest = true,
    checkdocs = :exports,
    linkcheck = true,
    html_disable_git = true, # disable source and edit links to github
    html_canonical = "https://juliangehring.github.io/MultipleTesting.jl/stable/",
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

deploydocs(
    repo = "github.com/juliangehring/MultipleTesting.jl.git",
    julia = "0.7",
    target = "build",
    deps   = nothing,
    make   = nothing
)
