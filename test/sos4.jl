using MomentTools, JuMP
using DynamicPolynomials
using MosekTools; opt = Mosek.Optimizer

set_optimizer(opt)

X = @polyvar x y 

f = 1-y+BigInt(57)//10
H = [x]
G = [1-x^2-y^2, x+y-1]

d = 2

s, P, Q, v, M = sos_decompose(f,H,G,X,d)

s1, P1, Q1, v1, M1 = esos_decompose(f, H, G, X, d)

wsos_decompose(s1)
