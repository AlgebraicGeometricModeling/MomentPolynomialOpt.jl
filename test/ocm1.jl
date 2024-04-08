using MomentPolynomialOpt, DynamicPolynomials, JuMP
using MosekTools; mpo_optimizer(Mosek.Optimizer, "QUIET"=> true)

X = @polyvar y u

d = 4

M = MOM.Model()

gamma = MOM.moment_variables(M, X, 2*d)
xi    = MOM.moment_variables(M, X, 2*d)

k = y
f = y*u

y0 = 0.5

g1 = 1 - y^2
g2 = 1 - u^2

MOM.add_constraint_nneg(M, g1, gamma)

MOM.add_constraint_nneg(M, g2, gamma)
MOM.add_constraint_unitmass(M, gamma)

MOM.add_constraint_nneg(M, g1, xi)
MOM.add_constraint_nneg(M, g2, xi)

N = 4
for i in 1:N
    @constraint(M, mmt(gamma, f^i)- mmt(gamma, y^i) ==0)
    @constraint(M, mmt(gamma, y) - mmt(gamma, y0) -   mmt(xi, f^i)+ mmt(xi, y^i) ==0)
end


@objective(M, Min, mmt(gamma, y))

optimize!(M)

s = get_series(M)

w, Gamma = decompose(s[1]) #gamma
v, Xi    = decompose(s[2]) #xi
