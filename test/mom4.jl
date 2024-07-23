using JuMP
using MomentPolynomialOpt
using DynamicPolynomials, MultivariateSeries

using MosekTools; mpo_optimizer(Mosek.Optimizer, "QUIET"=>true)

#using CSDP; mpo_optimizer(CSDP.Optimizer)

X = @polyvar x y

p =  x^2+2*x*y+1

d = 4

M = MOM.Model()

mu = moments(M, X, 2*d, :PRB)

@constraint(M, dot(mu, p) == 0)

MOM.set_objective_ncl(M, mu)

JuMP.optimize!(M)

s     = get_series(M)
w, Xi = get_measure(M)
