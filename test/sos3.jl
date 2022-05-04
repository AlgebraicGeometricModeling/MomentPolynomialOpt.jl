using MomentTools, JuMP
using DynamicPolynomials
using MosekTools; opt = Mosek.Optimizer
set_optimizer(opt)

X = @polyvar x y z u v

f = x
H = [x - u - v , y - u^2+u*v-v^2-1, z - u^3-v^3]
G = [u - u^2, v - v^2]

d = 2

M0 = SOS.Model(:sup,f,H,G,X,d, opt)
v0, M0 = optimize(M0)


M1 = SOS.Model(:sup,f,H,G,X,d)
v1, M1 = optimize(M1)

