using JuMP
using MomentTools
using DynamicPolynomials, MultivariateSeries

using MosekTools
optimizer = Mosek.Optimizer

#using CSDP; optimizer = CSDP.Optimizer

set_optimizer(optimizer)

X = @polyvar x y

p =  x^2+2*x*y+1


d = 4
M = MOM.Model(X, d; nu=2)

s0 = MultivariateSeries.dual(p)
L  = monomials(X,0:2)

MOM.constraint_moments(M, [(m=>s0[m]) for m in L])
MOM.set_objective_tv(M)

v, M = optimize(M)

s = get_series(M)

w, Xi = get_measure(M)
