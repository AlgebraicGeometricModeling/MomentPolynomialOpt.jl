using MomentTools
using DynamicPolynomials
using JuMP

using MosekTools
optimizer = Mosek.Optimizer
mmt_optimizer(optimizer)

X = @polyvar x y

p1 = x*y-1
p2 = x^2-2
g = 2*x-x^2

f = y

d = 3

M = SOS.Model(:inf, f, [p1,p2], [g], X, d)
v, M = optimize(M)
