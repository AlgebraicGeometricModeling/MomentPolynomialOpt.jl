module MomentTools

using DynamicPolynomials
using MultivariateSeries
using JuMP
using Dualization
using LinearAlgebra
using Combinatorics

MMT = Dict{Symbol,Any}( :dual => true)
export MMT

include("mommodel.jl")
#include("constraints.jl")
#include("objective.jl")

include("optimize.jl")
include("polar.jl")
include("annihilator.jl")

include("sosmodel.jl")
include("exact_decompose.jl")


end
