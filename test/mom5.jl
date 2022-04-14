using JuMP
using MomentTools
using DynamicPolynomials
using MultivariateSeries
using LinearAlgebra
using Dualization

# using Plots
# function plotmeas(w,Xi)
#     plot(real(Xi[1,:]), real(Xi[2,:]), seriestype = :scatter, zcolor = abs.(w), m = (:heat, 0.8, Plots.stroke(1, :black)))
# end

using MosekTools
if haskey(ENV,"QUIET")
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);
else
    optimizer = Mosek.Optimizer
end

#using CSDP; optimizer = CSDP.Optimizer


#MOM.set_optimizer(Dualization.dual_optimizer(optimizer))
MOM.set_optimizer(optimizer)
X = @polyvar x y

lebesgue(i,j) = ((1-(-1)^(i+1))/(i+1))*((1-(-1)^(j+1))/(j+1))

d = 10

M = MOM.Model(X, d; nu=2)

p1 = 1-x^2-y^2
# p1 * mu_1 >= 0
MOM.constraint_nneg(M, 1, 1-x^2-y^2)

q1 = 1-x^2
q2 = 1-y^2
# q1 * mu_2 >= 0, q2 * mu_2 >=0
MOM.constraint_nneg(M, 2, 1-x^2, 1-y^2 )

# monomials of degree <= 2*d
L = monomials(X, seq(0:2*d))

# <1*mu_1, m> + <1*mu_2, m> = leb_mom(m)
MOM.constraint_moments(M,
                   [(m=>lebesgue(exponents(m)...)) for m in L],
                   [1,1] )

# sup  <1*mu_1,1>  
MOM.set_objective(M, "sup", 1, 1)

optimize!(M)

v = objective_value(M)

println("Approximate volume: ", v)

# s = get_series(M)
# w, Xi = get_measure(M)

# sp = p1*s[1]
# L1 = monomials(X, seq(0:d-1)); H1 =  hankel(sp, L1, L1); vp1 = eigvals(H1)
# L2 = monomials(X, seq(0:d)); H2 =  hankel(s[1], L2, L2); vp2 = eigvals(H2)

#w, Xi = decompose(get_series(M)[1])
#plotmeas(w,Xi)

