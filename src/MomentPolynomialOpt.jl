module MomentPolynomialOpt

using DynamicPolynomials
#using MultivariateSeries
using JuMP
using Dualization
using LinearAlgebra
using Combinatorics

MPO = Dict{Symbol,Any}( :optimizer => nothing )
export MPO


include("moments.jl")

include("MOM/Model.jl")

include("SOS/Model.jl")
include("SOS/MaxEigenModel.jl")
include("SOS/MinEllipsoidModel.jl")

include("optimize.jl")
include("minimizers.jl")
include("polar.jl")
include("annihilator.jl")

include("sos_decompose.jl")

end



