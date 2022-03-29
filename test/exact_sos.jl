using LinearAlgebra, MomentTools, DynamicPolynomials, JuMP
#using MosekTools; optimizer=Mosek.Optimizer
using CSDP; MMT[:optimizer] = CSDP.Optimizer

#include("../src/exact_decompose.jl")


X = @polyvar x

f = (x^3-2)

g = x-124//100
#g = x^2 -x^5/4 + f*x^2/4

Q,q = exact_decompose(g,f,rounding=2, verbose=true)
