using MomentTools
using DynamicPolynomials
using JuMP

X = @polyvar x y

p1 = x*y-1
p2 = x^2-2
q  = 2*x-x^2
using MosekTools;optimizer = Mosek.Optimizer
#using MosekTools; optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);

using CSDP; optimizer = CSDP.Optimizer

d = 3

m = MOM.Model(X, d, optimizer)

MOM.add_constraint_measure(m)
MOM.add_constraint_zero(m, p1, p2)
MOM.add_constraint_nneg(m,q)
MOM.objective(m, y)

optimize!(m)

v  = objective_value(m)
y  = value(m)
L  = m[:monomials]
s  = getseries(m)
Xi = getminimizers(m)
