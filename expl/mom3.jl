using MomentTools
using DynamicPolynomials
using MultivariateSeries

using MosekTools;optimizer = Mosek.Optimizer
#using MosekTools; optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);

#using CSDP; optimizer = CSDP.Optimizer

X  = @polyvar x y
q1 = 1-x^2-y^2
q2 = x^3-y^2

v, m = minimize(x, [], [q1,q2], X, 4, optimizer)

w, Xi = getmeasure(m)



