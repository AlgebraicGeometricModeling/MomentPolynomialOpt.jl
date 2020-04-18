module MomentTools

using JuMP
using MultivariateSeries

include("mommodel.jl")
include("mosmodel.jl")
include("optimize.jl")

export MOM, getseries, getminimizers

end
