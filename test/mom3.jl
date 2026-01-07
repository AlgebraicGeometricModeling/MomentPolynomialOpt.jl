using MomentPolynomialOpt, DynamicPolynomials

using JuMP, MosekTools

mpo_optimizer(Mosek.Optimizer, "QUIET" => true)

X  = @polyvar x y
q1 = 1-x^2-y^2
q2 = x^3-y^2

pop = [(y,"inf"), (q1,">=0"), (q2, ">=0")]

v, M = optimize(pop, X, 4)

w, Xi = get_measure(M)



