using DynamicPolynomials
using MosekTools;optimizer = Mosek.Optimizer
using MomentTools

#using MosekTools; optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);
#using CSDP; optimizer = CSDP.Optimizer

X  = @polyvar x1 x2

e1 = x1^2-2
e2 = (x2^2-3)*(x1*x2-2)

p1 = x1
p2 = 2-x2

v, m = minimize(-x1, [e1, e2], [p1,p2], X, 3, optimizer)

println("value: ",v)

Xi = getminimizers(m)
