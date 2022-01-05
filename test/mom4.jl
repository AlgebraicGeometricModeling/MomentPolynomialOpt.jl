using JuMP
using MomentTools
using DynamicPolynomials
using MultivariateSeries
using MosekTools
if haskey(ENV,"QUIET")
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);
else
    optimizer = Mosek.Optimizer
end 
#using CSDP; optimizer = CSDP.Optimizer

X = @polyvar x y

p =  x^2+2*x*y+1



d = 4
M = MOM.Model(X, d; nu=2)
set_optimizer(M,optimizer)

s0 = MultivariateSeries.dual(p)
L = monomials(X,seq(0:2))

constraint_moments(M, [(m=>s0[m]) for m in L])

set_objective_tv(M)

v = optimize(M)[1]

s = get_series(M)

w, Xi = get_measure(M)
