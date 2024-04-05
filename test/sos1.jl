using MomentPolynomialOpt
using DynamicPolynomials
using JuMP

using MosekTools; mpo_optimizer(Mosek.Optimizer, "QUIET" =>true)


X = @polyvar x y

p1 = x*y-1
p2 = x^2-2
g = 2*x-x^2

f = y

d = 3

M = SOS.Model(:inf, f, [p1,p2], [g], X, d)
v, M = optimize(M)
