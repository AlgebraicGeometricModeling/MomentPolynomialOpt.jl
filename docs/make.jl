using Documenter
using MomentPolynomialOpt, JuMP

dir = "mrkd"
Expl = map(file -> joinpath("expl", file),
           filter(x ->endswith(x, ".md"), readdir(dir*"/expl")))
Code = map(file -> joinpath("code", file),
           filter(x ->endswith(x, ".md"), readdir(dir*"/code")))

makedocs(
    sitename = "MomentPolynomialOpt",
    authors = "L. Baldi, B. Mourrain",
    modules = [MomentPolynomialOpt],
    build = "MomentPolynomialOpt.jl/docs",
    source = "mrkd",
    pages = Any[
        "Home" => "index.md",
        "Tutorials" => Expl,
        "Manual" => Code,
        "About the package"  => "package.md",
    ],
    repo = Remotes.GitHub("AlgebraicGeometricModeling", "MomentPolynomialOpt.jl"),
    doctest = false,
)

#deploydocs(
#    deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
#    repo = "github.com/AlgebraicGeometricModeling/MomentPolynomialOpt.jl.git"
#)
