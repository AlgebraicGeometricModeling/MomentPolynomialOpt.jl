using MomentTools
using DynamicPolynomials
using JuMP

X = @polyvar x y

p1 = x*y-1
p2 = x^2-2

using MosekTools;optimizer = Mosek.Optimizer
#using MosekTools; optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);

#using CSDP; optimizer = CSDP.Optimizer

d = 4

m = MOM.Model(X, d, optimizer)

MOM.add_constraint_zero(m, p1, p2)
MOM.add_constraint_nneg(m, 4-x^2)
MOM.objective(m, x)

optimize!(m.model)
v = objective_value(m.model)

s = getseries(m)

Xi = getminimizers(m)
