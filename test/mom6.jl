using JuMP
using Combinatorics
using MomentTools
using DynamicPolynomials
using MultivariateSeries
using LinearAlgebra

using MosekTools
if haskey(ENV,"QUIET")
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true);
else
    optimizer = Mosek.Optimizer
end 

X = @polyvar x y

# Minimize on a compact singular semialgebraic set with polar_minimize
h = [] #zero constraints
g = [x^3-y^2, 1-x^2-y^2] #nonneg constraints
f = x #objective function
d = 10 #degree of approximation

v, m = polar_minimize(f, h, g, X, d, optimizer)
println(v)

# Minimize on a noncompact singular semialgebraic set with polar_minimize
h = [] #zero constraints
g = [x^3-y^2] #nonneg constraints
f = x #objective function
d = 10 #degree of approximation

v, m = polar_minimize(f, h, g, X, d, optimizer)
println(v)

# Minimize on a compact singular semialgebraic set with polar_minimize, adding the preordering
h = [] #zero constraints
g = [x^3-y^2, 1-x] #nonneg constraints
f = x #objective function
d = 10 #degree of approximation

pg = preordering(g)

v, m = polar_minimize(f, h, pg, X, d, optimizer)
println(v)

#Compare without addind the product: if it is possible, better to avoid to add it
w, n = polar_minimize(f, h, g, X, d, optimizer)
println(w)
println(v)
