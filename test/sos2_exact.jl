using LinearAlgebra, MomentPolynomialOpt, DynamicPolynomials, JuMP

#=
using MosekTools; opt=Mosek.Optimizer
if haskey(ENV,"QUIET")
    opt = optimizer_with_attributes(opt, "QUIET" => true);
end
=#
#using CSDP; opt=CSDP.Optimizer; mpo_optimizer(opt)
using MosekTools; mpo_optimizer(Mosek.Optimizer, "QUIET"=>true)

X = @polyvar x

h = x^3-2

f = x+0
#g = x^2 -x^5/4 + f*x^2/4

WS, P, v, M = sos_decompose(f, [h], [], X, 3; round=1, exact=true)

