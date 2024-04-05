using DynamicPolynomials, MomentPolynomialOpt, JuMP

using MosekTools; mpo_optimizer(Mosek.Optimizer,  "QUIET" => true)


X = @polyvar x y z

f = 1
H = [-y^2*z + x^2]
G = [-x^2 - y^2 - z^2 - 4*z - 3]

d = 3

v, M = optimize(:inf, f, H, G, X, d)

S = get_series(M)

K, I, P, B = annihilator(S,2)
