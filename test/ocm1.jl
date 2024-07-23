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

constraint_nneg(M, gamma, g1)
constraint_nneg(M, gamma, g2)

constraint_nneg(M, xi, g1)
constraint_nneg(M, xi, g2)

N = 4
for i in 1:N 
    @constraint(M, dot(gamma, f^i)- dot(gamma, y^i) ==0)
    @constraint(M, dot(gamma, y^i) - dot(gamma, y0^i) - dot(xi, f^i) + dot(xi, y^i) ==0)
end


@objective(M, Min, dot(gamma, y))

optimize!(M)

s = get_series(M)

v, Xi    = decompose(s[2]) #xi
w, Gamma = decompose(s[1]) #gamma
