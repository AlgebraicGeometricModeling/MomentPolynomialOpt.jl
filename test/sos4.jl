using MomentPolynomialOpt, JuMP
using DynamicPolynomials
#using MosekTools; opt = Mosek.Optimizer
using CSDP; opt = CSDP.Optimizer

mmt_optimizer(opt)

X = @polyvar x y 

f = 1-y+57//10
H = [x]
G = [1-x^2-y^2, x+y-1]

d = 2

WS0, WS, P, v, M = sos_decompose(f,G,H,X,d)

#s1, P1, Q1, v1, M1 = esos_decompose(f, H, G, X, d)

res = f - wsos(WS0) - dot(wsos.(WS), G) - dot(P, H)
println(" ==== ", norm(coefficients(res)))

X1 = @polyvar x


f = x^2+x+1
g = x
h = x^4-1

WS1, WS11, P,v, M = sos_decompose(f,[g],[h], X1, 2)
res = f - wsos(WS1) - wsos.(WS11)[1]*g - h*P[1]
println(" ==== ", norm(coefficients(res)))
