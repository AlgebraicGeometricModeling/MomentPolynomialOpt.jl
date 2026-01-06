export MomentTensorDecomposition

module MomentTensorDecomposition

# Export the main function to use in examples
export symm_tens_decomp

# List necessary dependencies
using MomentPolynomialOpt, DynamicPolynomials, MultivariateSeries, LinearAlgebra
using JuMP
using TensorDec

# Private functions starting with an underscore

function _symm_tens_decomp(X, l, F0, rescaling, use_kernel, ::Val{:real})
    println("Performing decomposition with real weights...")
    Y = vec(X[2:end])
    d = maxdegree(F0)
    F0 = F0(vcat(X[1], vec(X[2:end]) ./ rescaling))
    s0 = series(F0, Y) 

    ## Define model ###
    M = MOM.Model()
    
    ### Define variables ###
    mu = moments(M, Y, 2*l, :PSD)
    mu_ = moments(M, Y, 2*l, :PSD)
    ### Add non-negativity constraints ###
    g1 = 1 - sum(Y.^2)
    MOM.add_constraint_nneg(M, mu, g1)
    MOM.add_constraint_nneg(M, mu_, g1)

     ## Add kernel constraints ###
    if use_kernel == true
        null_polys = []
        for i in 1:Int(floor(d/2))
            L1 = monomials(Y, 0:i)
            L2 = monomials(Y, 0:i)
            H = hankel(s0, L1, L2)
            append!(null_polys, nullspace(H)'*vec(L2))
        end
        for p in null_polys
            MOM.add_constraint_zero(M, mu, p)
            MOM.add_constraint_zero(M, mu_, p)
        end
    end

     ### Add moment constraints ### 
    for (m,c) in zip(monomials(s0), coefficients(s0))
        @constraint(M, mmt(mu, m) - mmt(mu_, m)== c)
    end

    ### Define objective ###
    k = sum(monomials(Y, 0:l).^2) 
    @objective(M, Min, mmt(mu, k) + mmt(mu_, k))

    ### Solve ###   
    optimize!(M);
    S1, S2 = get_series(M);
    ### Decompose and compare generated tensor to the original ###
    w1, X1 = decompose(S1) 
    w2, X2 = decompose(S2) 

    if isempty(X1) && isempty(X2)
        error("Decomposition failed. Try modifying the inputs: e.g., changing the rescaling.")
    end

    X1 = vcat(fill(1.,size(X1,2))', X1)
    X2 = vcat(fill(1.,size(X2,2))', X2)
    F10 = MultivariateSeries.dual(s0, d)

    if !isempty(X2) && !isempty(X1)
        F1 = tensor(w1, X1, vcat(1, Y), d) + tensor(-w2, X2, vcat(1, Y), d) 
    elseif !isempty(X1)
        F1 = tensor(w1, X1, vcat(1, Y), d)
    elseif !isempty(X2)
        F1 = -tensor(w2, X2, vcat(1, Y), d)
    end
    norm = norm_apolar(F10 - F1)
    L = length(w1)+length(w2)

    ## Scale back the y and z variables ##
    X1[2:length(Y),:] = X1[2:length(Y),:].*rescaling 
    X2[2:length(Y),:] = X2[2:length(Y),:].*rescaling    

    if norm > 1
        println("Apolar distance to input polynomial is high. Try changing the input parameters.")
    end

    Xi = hcat(X1, X2)
    w = vcat(w1, -w2)

    return norm, L, Xi, w
end

function _symm_tens_decomp(X, l, F0, rescaling, use_kernel, ::Val{:positive})
    println("Performing decomposition with positive weights...")
    Y = vec(X[2:end])
    d = maxdegree(F0)
    F0 = F0(variables(F0) => vcat(X[1], Y./rescaling))
    s0 = series(F0, Y) 

    ### Define model ###
    M = MOM.Model()
    
    ### Define variables ###
    mu = moments(M, Y, 2*l, :PSD)
    ### Add non-negativity constraints ###
    g1 = 1 - sum(Y.^2)
    MOM.add_constraint_nneg(M, mu, g1)

    ### Add kernel constraints ###
    if use_kernel == true
        null_polys = []
        for i in 1:Int(floor(d/2))
            L1 = monomials(Y, 0:i)
            L2 = monomials(Y, 0:i)
            H = hankel(s0, L1, L2)
            append!(null_polys, nullspace(H)'*vec(L2))
        end
        for p in null_polys
            MOM.add_constraint_zero(M, mu, p)
        end
    end

    ### Add moment constraints ### 
    for (m,c) in zip(monomials(s0), coefficients(s0))
        @constraint(M, mmt(mu, m)== c)
    end

    ### Define objective ###
    k = sum(monomials(Y, 0:l).^2) 
    @objective(M, Min, mmt(mu, k))

    optimize!(M);
    S = get_series(M)
    w, Xi = decompose(S) 
    if isempty(Xi)
        error("Decomposition failed. Try modifying the inputs, or using a real decomposition.")
    end
    Xi = vcat(fill(1.,size(Xi,2))', Xi)
    F10 = MultivariateSeries.dual(s0, d)
    F1 = tensor(w, Xi, vcat(1, Y), d)
    norm = norm_apolar(F10 - F1)
    L = length(w)

    ## Scale back the y and z variables ##
    Xi[2:length(Y),:] = Xi[2:length(Y),:].*rescaling

    if norm > 1
        println("Apolar distance to input polynomial is high. Try changing the input parameters.")
    end



    return norm, L, Xi, w
end

# --- Public function ---

"""
Performs symmetric tensor decomposition of the polynomial F0.
Dispatches to the correct operation based on the 'weight_type'.
"""
function symm_tens_decomp(X, l, F0; rescaling=1, use_kernel=true, weight_type=:real)
    return _symm_tens_decomp(X, l, F0, rescaling, use_kernel, Val(weight_type))
end

end