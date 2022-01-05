ENV["QUIET"]= true


dir=pwd()

F0 = filter(x -> endswith(x, ".jl") && !startswith(x,"runtests"), readdir(dir))

for f in F0
    try
        t = @elapsed include(f)
        @info "\033[96m$f\033[0m   $t(s)"
    catch
        @warn "problem with $f"
    end
end

# using NBInclude
# F1 = filter(x ->endswith(x, ".ipynb"), readdir(dir*"/../expl"))
# for f in F1
#     @nbinclude(joinpath(dir*"/../expl", f))
# end
