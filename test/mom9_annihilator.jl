using DynamicPolynomials, MomentPolynomialOpt, PMO
using JuMP, MosekTools

import DynamicPolynomials: maxdegree
function DynamicPolynomials.maxdegree(i::Int64) return 0 end

mmt_optimizer(Mosek.Optimizer)


X = @polyvar x y z

f = 1
eq= [-y^2*z + x^2]
nn= [-x^2 - y^2 - z^2 - 4*z - 3]

d = 3

v, M = optimize(:inf, f, eq, nn, X, d)

S = get_series(M)[1]

K, I, P, B = annihilator(S,2)
