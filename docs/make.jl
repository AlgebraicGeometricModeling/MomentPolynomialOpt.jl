using Documenter
using MomentTools, JuMP

dir = "mrkd"
Expl = map(file -> joinpath("expl", file),
           filter(x ->endswith(x, ".md"), readdir(dir*"/expl")))
Code = map(file -> joinpath("code", file), filter(x ->endswith(x, ".md"), readdir(dir*"/code")))

makedocs(
    sitename = "MomentTools",
#    format = Documenter.HTML(prettyurls = false),
    authors = "L. Baldi, B. Mourrain",
    modules = [MomentTools],
    build = "MomentTools.jl/docs",
    source = "mrkd",
    pages = Any[
        "Home" => "index.md",
        "Functions & types" => Code,
        "Example" => Expl,
        "About the package"  => "package.md",
    ],
    repo = "https://github.com/AlgebraicGeometricModeling/MomentTools.jl",
    doctest = false,
)

deploydocs(
    deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo = "github.com/AlgebraicGeometricModeling/MomentTools.jl.git"
)
