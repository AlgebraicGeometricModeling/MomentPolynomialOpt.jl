using MomentPolynomialOpt, DynamicPolynomials, MultivariateSeries, JuMP
using MosekTools; mpo_optimizer(Mosek.Optimizer, "QUIET"=> true)

X = @polyvar y u

d = 4

M = MOM.Model()

gamma = moments(M, X, 2*d, :PRB)
xi    = moments(M, X, 2*d, :PSD)

k = y
f = y*u

y0 = 0.5

g1 = 1 - y^2
g2 = 1 - u^2

MOM.add_constraint_nneg(M, g1, gamma)
MOM.add_constraint_nneg(M, g2, gamma)


MOM.add_constraint_nneg(M, g1, xi)
MOM.add_constraint_nneg(M, g2, xi)

N = 4
for i in 1:N
    @constraint(M, mmt(gamma, f^i)- mmt(gamma, y^i) ==0)
    @constraint(M, mmt(gamma, y^i) - mmt(gamma, y0^i) -   mmt(xi, f^i)+ mmt(xi, y^i) ==0)
end


@objective(M, Min, mmt(gamma, y))

optimize!(M)

s = get_series(M)

v, Xi    = decompose(s[2]) #xi
w, Gamma = decompose(s[1]) #gamma
