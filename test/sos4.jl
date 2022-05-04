using MomentTools, JuMP
using DynamicPolynomials
using MosekTools; opt = Mosek.Optimizer

set_optimizer(opt)

X = @polyvar x y 

f = 1-y+0.1
H = [x]
G = [1-x^2-y^2, x+y-1]

d = 2

s, p, q, v, M = sos_decompose(f,H,G,X,d)



