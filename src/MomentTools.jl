module MomentTools

using JuMP
using MultivariateSeries

include("mommodel.jl")
include("optimize.jl")

export MOM, getseries, getminimizers, getmeasure

end
