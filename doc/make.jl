using Documenter
using MomentTools

dir = "mrkd"
Expl = map(file -> joinpath("expl", file),
           filter(x ->endswith(x, "md"), readdir(dir*"/expl")))
#Code = map(file -> joinpath("code", file), filter(x ->endswith(x, "md"), readdir(dir*"code")))

makedocs(
    sitename = "MomentTools",
    format = Documenter.HTML(),
    authors = "B. Mourrain",
    modules = [MomentTools],
    build = "html",
    source = dir,
    pages = Any[
        "Home" => "index.md",
        "Example" => Expl,
#        "Functions & types" => Code
    ],
    doctest = false
)

deploydocs(
    repo = "gitlab.inria.fr/AlgebraicGeometricModeling/MomentTools.jl.git"
)
