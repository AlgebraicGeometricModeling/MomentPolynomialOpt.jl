module MomentTools

using DynamicPolynomials
using MultivariateSeries
using JuMP
using Dualization
using LinearAlgebra
using Combinatorics

include("mommodel.jl")
include("constraints.jl")
include("objective.jl")
include("optimize.jl")
include("polar.jl")
include("annihilator.jl")
include("exact_decompose.jl")

MMT = Dict{Symbol,Any}()
export MMT

end
