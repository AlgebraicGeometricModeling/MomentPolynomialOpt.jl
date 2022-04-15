using MomentTools
using DynamicPolynomials
using JuMP

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
q  = 2*x-x^2
f  = y

d = 3

M = MOM.Model(X, d, optimizer)
#mom_optimizer(M, optimizer)

MOM.constraint_unitmass(M)
MOM.constraint_zero(M, p1, p2)
MOM.constraint_nneg(M, q)
MOM.set_objective(M, "inf", f)

v = optimize(M)[1]

s     = get_series(M)
w, Xi = get_measure(M)
