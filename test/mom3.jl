using MomentTools
using DynamicPolynomials
using MultivariateSeries
using JuMP, MosekTools

#if haskey(ENV,"QUIET")
#    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);
#else
optimizer = Mosek.Optimizer
#end 
set_optimizer(optimizer)

X  = @polyvar x y
q1 = 1-x^2-y^2
q2 = x^3-y^2

pop = [(x,"inf"), (q1,">=0"), (q2, ">=0")]

v, M = MOM.optimize(pop, X, 4)

w, Xi = get_measure(M)



