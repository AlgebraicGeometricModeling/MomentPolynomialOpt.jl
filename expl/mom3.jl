using MomentTools
using DynamicPolynomials
using MultivariateSeries

using MosekTools;optimizer = Mosek.Optimizer
#using MosekTools; optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);

#using CSDP; optimizer = CSDP.Optimizer

X  = @polyvar x y
q1 = 1-x^2-y^2
q2 = x^3-y^2

pop = [(x,"inf"), (q1,">=0"), (q2, ">=0")]
v, m = optimize(pop, X, 4, optimizer)

w, Xi = getmeasure(m)



