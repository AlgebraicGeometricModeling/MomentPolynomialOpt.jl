using DynamicPolynomials, MomentTools, PMO
using JuMP; MosekTools

import DynamicPolynomials: maxdegree

function DynamicPolynomials.maxdegree(i::Int64) return 0 end

L = PMO.table()

for i in 1:length(L)
    local p = L[i]
    if p[:type] == "polynomial" && p[:nvar] <= 6 
        try
            print("---> ",p[:name], ": ")
            local P = vec(p)
            local d = max([maxdegree(P[i][1]) for i in 1:length(P)]...)
            
            local t = @elapsed v, M = optimize(P, variables(p), div(d+1,2),
                        optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true))

            println(" nv: ", p[:nvar],"  d: ", d, "  o: ", div(d+1,2), "  v*: ", v, "  ",t,"(s)")
        
        catch
            @warn "error while building the optimization model"
        end
    end
end
