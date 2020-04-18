using MomentTools
using DynamicPolynomials

using MosekTools;optimizer = Mosek.Optimizer
#using MosekTools; optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);

#using CSDP; optimizer = CSDP.Optimizer

X  = @polyvar x y
p1 = 1-x^2-y^2
p2 = x^3-y^2

v, m = minimize(x, [], [p1,p2], X, 4, optimizer)

getminimizers(m)

