using MomentTools
using DynamicPolynomials
using JuMP
#using Dualization
using MosekTools
if haskey(ENV,"QUIET")
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);
else
    optimizer = Mosek.Optimizer
end
set_optimizer(optimizer)

X = @polyvar x y

p1 = x*y-1
p2 = x^2-2
g = 2*x-x^2

f = y

d = 3

v, M = SOS.optimize(:Min, f, [p1,p2], [g], X, d)
