using MomentPolynomialOpt, DynamicPolynomials, TensorDec
using MosekTools; mpo_optimizer(Mosek.Optimizer, "QUIET" =>true)
#using CSDP; mpo_optimizer(CSDP.Optimizer)

println("--- Positive Weight Decomposition Example (tensor2.jl) ---")

## Input parameters ##
X = @polyvar x y z
l = 6
rescaling = 20
weight_type = :positive

# Define the polynomial to be decomposed
F0 = -1549440*x*y*z^3 + 2417040*x*y^2*z^2 + 166320*x^2*y*z^2 - 829440*x*y^3*z - 5760*x^3*y*z - 222480*x^2*y^2*z +
    38*x^5 - 497664*y^5 - 1107804*z^5 - 120*x^4*y + 180*x^4*z + 12720*x^3*y^2 + 8220*x^3*z^2 - 34560*x^2*y^3 - 59160*x^2*z^3 +
    831840*x*y^4 + 442590*x*z^4 - 5591520*y^4*z + 7983360*y^3*z^2 - 9653040*y^2*z^3 + 5116680*y*z^4
F0 = differentiate(F0, x) / 5

# Call function from the Tensor Decomposition module
Xi, w = MomentTensorDecomposition.symm_tens_decomp(X, l, F0, rescaling = rescaling, weight_type = weight_type)

println("\n--- Results ---")

## Construct the approximated tensor ##
F1 = tensor(w, Xi, vcat(1, vec(X[2:end])), maxdegree(F0))
## Homogenize it ##
Fh = sum(coefficient(F1,m)*m*x^(maxdegree(F0)-maxdegree(m)) for m in monomials(F1))
## Compute apolar norm between F0 and Fh ##
norm = norm_apolar(F0 - Fh)

## Print results ##
println("Apolar norm: ", norm)
println("Decomposition length: ", length(w))
println("Decomposition weights: ")
println(join([rpad(string(round(x, digits=3)), 10) for x in w], ""))
println("Decomposition points: ")
for row in eachrow(Xi)
    println(join([rpad(string(round(x, digits=3)), 10) for x in row], ""))
end
println("-----------------\n")
