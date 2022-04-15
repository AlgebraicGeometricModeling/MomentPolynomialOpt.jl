using JuMP, DynamicPolynomials, MomentTools

using MosekTools
if haskey(ENV,"QUIET")
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);
else
    optimizer = Mosek.Optimizer
end

#using CSDP; optimizer = CSDP.Optimizer
set_optimizer(optimizer)

X  = @polyvar x1 x2

e1 = x1^2-2
e2 = (x2^2-3)*(x1*x2-2)

p1 = x1
p2 = 2-x2

f = x1
v, M = MOM.minimize(f, [e1, e2], [p1, p2], X, 4)

Xi = get_minimizers(M)
