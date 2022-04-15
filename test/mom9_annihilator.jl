using DynamicPolynomials, MomentTools, PMO
using JuMP, MosekTools

import DynamicPolynomials: maxdegree
function DynamicPolynomials.maxdegree(i::Int64) return 0 end

set_optimizer(Mosek.Optimizer)


X = @polyvar x y z

f = 1
eq= [-y^2*z + x^2]
nn= [-x^2 - y^2 - z^2 - 4*z - 3]

#v, M = minimize(f, eq, nn, X, 2,
#                optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true))

d = 3

v, M = MOM.optimize(:inf, 1, eq, nn, X, d)

S = get_series(M)[1]

K, I, P, B = annihilator(S,2)
