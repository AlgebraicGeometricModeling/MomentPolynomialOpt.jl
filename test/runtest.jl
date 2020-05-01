using NBInclude
dir=pwd()
F0 = filter(x ->startswith(x, "mom"), readdir(dir))

for f in F0
    include(f)
end

F1 = filter(x ->endswith(x, ".ipynb"), readdir(dir*"/../expl"))

for f in F1
    @nbinclude(joinpath(dir*"/../expl", f))
end
