using JuMP
using MomentTools
using DynamicPolynomials
using MultivariateSeries

X = @polyvar x y

p =  x^2+2*x*y+1

using MosekTools;optimizer = Mosek.Optimizer
#using MosekTools; optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);

#using CSDP; optimizer = CSDP.Optimizer

d = 4
M = MOM.Model(X, d; nu=2)
set_optimizer(M,optimizer)

s = MultivariateSeries.dual(p)
L = monomials(X,seq(0:2))

constraint_moments(M, [(m=>s[m]) for m in L])

objective_tv(M)

optimize!(M)

v = objective_value(M)
s = getseries(M)

w, Xi = getmeasure(M)
