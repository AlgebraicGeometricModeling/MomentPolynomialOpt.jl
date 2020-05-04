using Documenter, DocumenterMarkdown
using MomentTools

dir = "mrkd"
Expl = map(file -> joinpath("expl", file),
           filter(x ->endswith(x, "md"), readdir(dir*"/expl")))
Code = map(file -> joinpath("code", file), filter(x ->endswith(x, "md"), readdir(dir*"code")))

makedocs(
    sitename = "MomentTools",
    format = Documenter.HTML(prettyurls = false),
    authors = "B. Mourrain",
    modules = [MomentTools],
    build = "html",
    source = dir,
    pages = Any[
        "Home" => "index.md",
        "Example" => Expl,
#        "Test" => ["a" => "try.html"],
#        "Functions & types" => Code
    ],
    doctest = false
)

deploydocs(
    deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo = "gitlab.inria.fr/AlgebraicGeometricModeling/MomentTools.jl.git"
)
