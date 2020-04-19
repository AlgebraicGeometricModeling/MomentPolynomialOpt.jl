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
m = MOM.Model(X, d, optimizer; nu=2)
MOM.add_constraint_moments(m, MultivariateSeries.dual(p))
MOM.objective_tv(m)

optimize!(m)

v = objective_value(m)
s = getseries(m)

w, Xi = getmeasure(m)
