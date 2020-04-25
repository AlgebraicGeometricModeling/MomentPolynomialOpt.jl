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

M = MOM.Model(X, d)
set_optimizer(M, optimizer)
constraint_unitmass(M)
constraint_zero(M, p1, p2)
constraint_nneg(M, q)
objective(M, y)

optimize!(M)

v  = objective_value(M)
y  = value(M)
L  = M[:monomials]

s  = getseries(M)
Xi = getminimizers(M)
