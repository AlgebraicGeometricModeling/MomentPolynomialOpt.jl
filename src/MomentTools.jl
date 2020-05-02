module MomentTools

using DynamicPolynomials
using MultivariateSeries
using JuMP
using Dualization

include("mommodel.jl")
include("constraints.jl")
include("objective.jl")
include("optimize.jl")


export MOM, getseries, getminimizers, getmeasure

end
