using MomentTools, DynamicPolynomials
using MosekTools, JuMP; opt = optimizer_with_attributes(Mosek.Optimizer, "QUIET"=> true)
mmt_optimizer(opt)

X = @polyvar x y

p1 = x*y-1
p2 = x^2-2
eq=[p1,p2]

q  = 2*x-x^2
ineq=[q]

f  = y

d = 3

#=
M = MOM.Model(X, d)
MOM.constraint_unitmass(M)
MOM.constraint_zero(M, p1, p2)
MOM.constraint_nneg(M, q)
MOM.set_objective(M, "inf", f)
MOM.dualize!(M)
=#

v, M  = optimize(:inf, f,eq, ineq, X, d)

s     = get_series(M)
w, Xi = get_measure(M)
