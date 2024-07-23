using JuMP
using MomentPolynomialOpt
using DynamicPolynomials
using MultivariateSeries
using LinearAlgebra

# using Plots
# function plotmeas(w,Xi)
#     plot(real(Xi[1,:]), real(Xi[2,:]), seriestype = :scatter, zcolor = abs.(w), m = (:heat, 0.8, Plots.stroke(1, :black)))
# end

using MosekTools; mpo_optimizer(Mosek.Optimizer, "QUIET" => true)


#using CSDP; optimizer = CSDP.Optimizer


X = @polyvar x y

function lebesgue(m)
    e = exponents(m); 
    return ((1-(-1)^(e[1]+1))/(e[1]+1))*((1-(-1)^(e[2]+1))/(e[2]+1))
end
d = 10

M = MOM.Model()


mu1 = moments(M,X,2*d,:PSD)

g1 = 1-x^2-y^2

# p1 * mu >= 0
constraint_nneg(M, mu1, g1)

mu2 = moments(M, X, 2*d, :PSD)

q1 = 1-x^2
q2 = 1-y^2

# q1 * mu_2 >= 0, q2 * mu_2 >=0
constraint_nneg(M, mu2, q1)
constraint_nneg(M, mu2, q2)

# monomials of degree <= 2*d
L = monomials(X, 0:2*d)

# <1*mu_1, m> + <1*mu_2, m> = leb_mom(m)

for m in L
    @constraint(M, dot(mu1, m) + dot(mu2,m) - lebesgue(m) == 0)
end

# sup  <1*mu_1,1>  
@objective(M, Max, dot(mu1,1) )

v, M = optimize(M)

#JuMP.optimize!(M)
#v = JuMP.objective_value(M)

println("Approximate volume: ", v)

S     = get_series(M)
#w, Xi = get_measure(M)


