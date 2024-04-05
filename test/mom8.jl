using MomentPolynomialOpt
using DynamicPolynomials
using JuMP

using MosekTools
optimizer = Mosek.Optimizer

mpo_optimizer(Mosek.Optimizer, "QUIET"=> true)

X = @polyvar x y

p1 = x*y-1
p2 = x^2-2
q  = 2*x-x^2

d = 3

v, M = optimize(:Inf, y, [p1,p2], [q], X, d)

s     = get_series(M)
w, Xi = get_measure(M)
