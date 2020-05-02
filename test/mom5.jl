using JuMP
using MomentTools
using DynamicPolynomials
using MultivariateSeries
using LinearAlgebra

# using Plots
# function plotmeas(w,Xi)
#     plot(real(Xi[1,:]), real(Xi[2,:]), seriestype = :scatter, zcolor = abs.(w), m = (:heat, 0.8, Plots.stroke(1, :black)))
# end

using MosekTools;optimizer = Mosek.Optimizer
#using MosekTools; optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);
#using CSDP; optimizer = CSDP.Optimizer

X = @polyvar x y

leb_mom(i,j) = ((1-(-1)^(i+1))/(i+1))*((1-(-1)^(j+1))/(j+1))

d = 10

M = MOM.Model(X, d, optimizer; nu=2)

p1 = 1-x^2-y^2
# p1 * mu_1 >= 0
constraint_nneg(M, 1, p1)

q1 = 1-x^2
q2 = 1-y^2
# q1 * mu_2 >= 0, q2 * mu_2 >=0
constraint_nneg(M, 2, q1, q2 )

# monomials of degree <= 2*d
L = monomials(X, seq(0:2*d))

# <1*mu_1, m> + <1*mu_2, m> = leb_mom(m)
constraint_moments(M,
                   [(m=>leb_mom(exponents(m)...)) for m in L],
                   collect(1:2), [1,1] )

# sup  <1*mu_1,1>  
objective(M, [1], [1], "sup")

v = optimize(M)[1]
println("Approximate volume: ", v)

# s = getseries(M)
# w, Xi = getmeasure(M)

# sp = p1*s[1]
# L1 = monomials(X, seq(0:d-1)); H1 =  hankel(sp, L1, L1); vp1 = eigvals(H1)
# L2 = monomials(X, seq(0:d)); H2 =  hankel(s[1], L2, L2); vp2 = eigvals(H2)

#w, Xi = decompose(getseries(M)[1])
#plotmeas(w,Xi)

