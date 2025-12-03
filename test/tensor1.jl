using MomentPolynomialOpt, DynamicPolynomials
using MosekTools; mpo_optimizer(Mosek.Optimizer, "QUIET" =>true)

println("--- Real Weight Decomposition Example (tensor1.jl) ---")

## Input parameters ##
X = @polyvar x y z t
l = 6
rescaling = 2
weight_type = :real

# Define the polynomial to be decomposed
F0 = -0.2489979301598193 * t^4 - 0.5714471586952264 * z * t^3 - 0.5050309495726817 * z^2 * t^2 + 0.20686591014734546 * z^3 * t +
    0.35848687448905797 * z^4 + 0.6290084096294283 * y * t^3 + 1.5863914541246662 * y * z * t^2 + 0.09102584099156069 * y * z^2 * t -
    0.7995202943157504 * y * z^3 - 1.000984626080145 * y^2 * t^2 - 0.7330971731383541 * y^2 * z * t + 0.48952697145340573 * y^2 * z^2 +
    0.3654972483140739 * y^3 * t + 0.04180994666887122 * y^3 * z - 0.06496525932745015 * y^4 + 0.8614765507794534 * x * t^3 +
    2.518149761933887 * x * z * t^2 + 0.11198414551091718 * x * z^2 * t - 1.2824172455614082 * x * z^3 - 3.198429373199247 * x * y * t^2 -
    2.280206948221652 * x * y * z * t + 1.51085249826687 * x * y * z^2 + 1.8431809395581764 * x * y^2 * t + 0.38190508335832307 * x * y^2 * z -
    0.6712764568310912 * x * y^3 - 2.2452396434195863 * x^2 * t^2 - 1.3839699114386361 * x^2 * z * t + 1.4125957148603887 * x^2 * z^2 +
    1.4670698710382422 * x^2 * y * t - 0.015953647647009017 * x^2 * y * z - 0.361365758387135 * x^2 * y^2 + 1.253344034777953 * x^3 * t +
    0.11846876725173883 * x^3 * z - 2.1207614308184404 * x^3 * y + 0.6141543193928954 * x^4

# Call the main function from the module
anorm, L, Xi, w = MomentTensorDecomposition.symm_tens_decomp(X, l, F0, rescaling = rescaling, weight_type = weight_type)

println("\n--- Results ---")
println("Apolar norm: ", anorm)
println("Decomposition length: ", L)
println("Decomposition weights: ")
println(join([rpad(string(round(x, digits=3)), 10) for x in w], ""))
println("Decomposition points: ")
for row in eachrow(Xi)
    println(join([rpad(string(round(x, digits=3)), 10) for x in row], ""))
end
println("-----------------\n")