using MomentTools
using DynamicPolynomials
using JuMP

using MosekTools
optimizer = Mosek.Optimizer

set_optimizer(optimizer)

X = @polyvar x y

p1 = x*y-1
p2 = x^2-2
q  = 2*x-x^2

d = 3

M = MOM.Model(:Inf, y, [p1,p2], [q], X, d)

#MOM.set_optimizer(M, optimizer)

v, M = optimize(M)

s     = get_series(M)
w, Xi = get_measure(M)
