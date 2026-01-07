using DynamicPolynomials, MomentPolynomialOpt

using MosekTools; mpo_optimizer(Mosek.Optimizer,  "QUIET" => true)
#using CSDP; mpo_optimizer(CSDP.Optimizer)

X = @polyvar x y

f = x + y + 3
G = [ y ]
H = [ x^2-1, y^2-x-2 ]

d = 2

L  = monomials(X, 0:d)

WS, P, v, M = sos_decompose(f, H, G, X, d; exact=true, round = 1 )

eq0 = f - wsos(WS[1]) - dot(G, wsos.(WS[2:end])) - dot(H,P)
println("reminder === ", eq0 )

Q0 = wsos(WS[1])
Q1 = wsos(WS[2])

