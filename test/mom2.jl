using JuMP, DynamicPolynomials, MomentTools
    
using MosekTools
optimizer = JuMP.optimizer_with_attributes(Mosek.Optimizer,"QUIET" => true)

#using CSDP; optimizer = CSDP.Optimizer
mmt_optimizer(optimizer)

X  = @polyvar x1 x2

e1 = x1^2-2
e2 = (x2^2-3)*(x1*x2-2)

p1 = x1
p2 = 2-x2

f = x1

v, M = minimize(f, [e1, e2], [p1, p2], X, 4)

Xi = get_minimizers(M)
