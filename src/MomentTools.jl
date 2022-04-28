module MomentTools

using DynamicPolynomials
using MultivariateSeries
using JuMP
using Dualization
using LinearAlgebra
using Combinatorics

MMT = Dict{Symbol,Any}( :optimizer => nothing )
export MMT

include("mom_model.jl")

include("minimizers.jl")
include("polar.jl")
include("annihilator.jl")

include("sos_model.jl")
include("sos_decompose.jl")
include("sos_exact.jl")


export optimize
#----------------------------------------------------------------------
"""
```julia
v, M = optimize(M)
```
Run the optimizer on the moment program `M` and output the objective_value `v` and the moment program `M`. If the optimization program has no value, it returns `nothing` and `M`.
"""
function optimize(M::JuMP.Model)
    JuMP.optimize!(M)
    if JuMP.has_values(M)
        return JuMP.objective_value(M), M
    else
        println("Solver status: ", JuMP.termination_status(M))
        return nothing, M
    end
end

import JuMP: set_optimizer

"""
```julia
set_optimizer(opt)
```
Define the default optimizer `opt` for the optimization problems created by MomentTools
"""
function JuMP.set_optimizer(opt)
    MMT[:optimizer] = opt 
end

end

