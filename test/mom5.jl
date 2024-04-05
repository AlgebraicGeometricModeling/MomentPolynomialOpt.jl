using JuMP
using MomentPolynomialOpt
using DynamicPolynomials
using MultivariateSeries
using LinearAlgebra

# using Plots
# function plotmeas(w,Xi)
#     plot(real(Xi[1,:]), real(Xi[2,:]), seriestype = :scatter, zcolor = abs.(w), m = (:heat, 0.8, Plots.stroke(1, :black)))
# end

using MosekTools
optimizer = Mosek.Optimizer

#using CSDP; optimizer = CSDP.Optimizer

mpo_optimizer(optimizer, "QUIET" => true)

X = @polyvar x y

lebesgue(i,j) = ((1-(-1)^(i+1))/(i+1))*((1-(-1)^(j+1))/(j+1))

d = 10

M = MOM.Model(X,d)


mu1 = M[:mu][1]

g1 = 1-x^2-y^2

# p1 * mu >= 0
MOM.add_constraint_nneg(M, g1, mu1)

mu2 = MOM.add_variable_moments(M, X, 2*d, :mu2)
MOM.add_constraint_nneg(M, mu2)


q1 = 1-x^2
q2 = 1-y^2
# q1 * mu_2 >= 0, q2 * mu_2 >=0
MOM.add_constraint_nneg(M, q1, mu2)
MOM.add_constraint_nneg(M, q2, mu2)

# monomials of degree <= 2*d
L = monomials(X, 0:2*d)

# <1*mu_1, m> + <1*mu_2, m> = leb_mom(m)

for m in L
    @constraint(M, mmt(mu1, m) + mmt(mu2,m) - lebesgue(exponents(m)...) == 0)
end

# sup  <1*mu_1,1>  
@objective(M, Max, mmt(mu1,1) )

JuMP.optimize!(M)
v = JuMP.objective_value(M)

println("Approximate volume: ", v)

s = get_series(M)
w, Xi = get_measure(M)


