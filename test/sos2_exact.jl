using LinearAlgebra, MomentTools, DynamicPolynomials, JuMP

#=
using MosekTools; opt=Mosek.Optimizer
if haskey(ENV,"QUIET")
    opt = optimizer_with_attributes(opt, "QUIET" => true);
end
=#
using CSDP; opt=CSDP.Optimizer

mmt_optimizer(opt)



#import MomentTools.exact_decompose
#include("../src/exact_decompose.jl")


X = @polyvar x

f = (x^3-2)

g = x - 0
#g = x^2 -x^5/4 + f*x^2/4

Q, q = esos_decompose(g,f,rounding=1, verbose=true)
