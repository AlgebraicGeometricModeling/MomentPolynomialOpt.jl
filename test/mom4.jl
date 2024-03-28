using JuMP
using MomentPolynomialOpt
using DynamicPolynomials, MultivariateSeries

using MosekTools; mmt_optimizer(Mosek.Optimizer, "QUIET"=>true)

#using CSDP; opt = CSDP.Optimizer



X = @polyvar x y

p =  x^2+2*x*y+1

d = 4

s0 = MultivariateSeries.dual(p)
L  = monomials(X,0:2)

M = MOM.Model(X, d; nu=2)
MOM.constraint_moments(M, [(m=>s0[m]) for m in L])
MOM.set_objective_tv(M)
MOM.dualize!(M,opt)

v, M = optimize(M)

s = get_series(M)

w, Xi = get_measure(M)
