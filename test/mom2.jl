using JuMP, DynamicPolynomials, MomentPolynomialOpt
    
using MosekTools

#using CSDP; mpo_optimizer(CSDP.Optimizer)
mpo_optimizer(Mosek.Optimizer,"QUIET" => true)

X  = @polyvar x1 x2

e1 = x1^2-2
e2 = (x2^2-3)*(x1*x2-2)

p1 = x1
p2 = 2-x2

f = x1

v, M = minimize(f, [e1, e2], [p1, p2], X, 4)

Xi = get_minimizers(M)
