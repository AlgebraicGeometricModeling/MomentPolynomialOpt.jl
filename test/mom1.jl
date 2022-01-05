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

X = @polyvar x y

p1 = x*y-1
p2 = x^2-2
q  = 2*x-x^2

d = 3

M = MOM.Model(X, d)
set_optimizer(M, optimizer)

constraint_unitmass(M)
constraint_zero(M, p1, p2)
constraint_nneg(M, q)
set_objective(M, "inf", y)

v = optimize(M)[1]

s     = get_series(M)
w, Xi = get_measure(M)
