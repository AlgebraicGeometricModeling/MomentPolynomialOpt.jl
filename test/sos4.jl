using MomentPolynomialOpt, JuMP
using DynamicPolynomials, LinearAlgebra
using MosekTools
mpo_optimizer(Mosek.Optimizer, "QUIET"=>true)
#using CSDP; opt = CSDP.Optimizer;mpo_optimizer(opt)

X = @polyvar x y 

f = 1-y+57//10
H = [x]
G = [1-x^2-y^2, x+y-1]

d = 2

WS0, P0, v, M = sos_decompose(f,H,G,X,d)


res0 = f - wsos(WS0[1]) - dot(wsos.(WS0[2:end]), G) - dot(P0, H)
println(" ==== ", norm(coefficients(res0)))


X1 = @polyvar x

f = x^2+x+1
G1 = [x]
H1 = [x^4-1]

v, M = minimize(f, H1, G1, X1, 3)
println("min:  ",v)

WS1, P1, v1, M = sos_decompose(f, H1, G1, X1, 2) #, round=1)


res1 = f - wsos(WS1[1]) - wsos(WS1[2])*G1[1] - H1[1]*P1[1]
println(" ==== ", res1)

