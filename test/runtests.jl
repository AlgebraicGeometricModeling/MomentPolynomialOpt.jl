ENV["QUIET"]= true


dir=pwd()

F0 = filter(x -> endswith(x, ".jl") && !startswith(x,"runtests"), readdir(dir))

E = String[]
T = Float64[]

for f in F0
    try
        s = @elapsed include(f)
        @info "\033[96m$f\033[0m   $s(s)"
        push!(T,s)
    catch
        @warn "problem with $f"
        push!(E,f)
    end
end

if length(E)>0
    @warn "problems with $E"
end

tt = sum(T)
@info "total time: $tt(s)"
