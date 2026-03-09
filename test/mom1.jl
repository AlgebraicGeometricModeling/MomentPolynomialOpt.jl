using MomentPolynomialOpt, DynamicPolynomials
#using MosekTools; mpo_optimizer(Mosek.Optimizer)
#using MosekTools; mpo_optimizer(Mosek.Optimizer, "QUIET"=> true)
using CSDP; mpo_optimizer(CSDP.Optimizer)

X = @polyvar x y

h1 = x*y-1
h2 = x^2-2
H=[h1,h2]

g = 2*x-x^2
G = [g]

f  = y

d = 3

v, M  = optimize(:inf, f, H, G, X, d)

s     = get_series(M)
w, Xi = get_measure(M)
