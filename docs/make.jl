using Documenter
using MultipleTesting

makedocs(modules = [MultipleTesting],
    sitename = "MultipleTesting.jl",
    format = Documenter.HTML(disable_git = true,  # disable source and edit links to github
        canonical = "https://juliangehring.github.io/MultipleTesting.jl/stable/"),
    doctest = true,
    checkdocs = :exports,
    linkcheck = true,
    pages = [
        "Home" => "index.md",
        "Library" => Any[
            "adjustment.md",
            "combination.md",
            "pi0.md",
            "higher-criticism.md",
            "models.md",
        ]
    ])

deploydocs(repo = "github.com/juliangehring/MultipleTesting.jl.git")