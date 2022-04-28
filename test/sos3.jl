using MomentTools, JuMP
using DynamicPolynomials
using MosekTools; opt = Mosek.Optimizer
set_optimizer(opt)

X = @polyvar x y z u v

f = x
H = [x - u - v , y - u^2+u*v-v^2-1, z - u^3-v^3]
G = [u - u^2, v - v^2]

d = 2

M = SOS.Model(:Inf,f,H,G,X,d, opt)


v0, M0 = optimize(M)


v1, M1 = SOS.optimize(:sup,f,H,G,X,d)

using CSDP; opt1= CSDP.Optimizer

v2, M2 = SOS.maximize(f,H,G,X,d, opt1)
