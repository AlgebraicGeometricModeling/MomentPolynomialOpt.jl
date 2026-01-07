export MomentTensorDecomposition

module MomentTensorDecomposition

# Export the main function to use in examples
export symm_tens_decomp

# List necessary dependencies
using MomentPolynomialOpt, DynamicPolynomials, MultivariateSeries, LinearAlgebra
using JuMP

#----------------------------------------------------------------------
#== 
```
_series_from_hpol(F, X0, Y, rescaling, d = maxdegree(F))
```
Convert a homogeneous polynomial `F` to a series representation by setting variable `X0` to 1,
rescaling variables `Y`, and normalizing coefficients by binomial factors.
==#
function _series_from_hpol(F, X0, Y, rescaling, d = maxdegree(F))
    P = subs(F,X0=>1, [y=>y/rescaling for y in Y]...)
    c = coefficients(P)
    m = monomials(P)
    for i in 1:length(c)
        e = exponents(m[i])
        c[i]/=binomial(d, e)
    end
    return MultivariateSeries.dual(P)
end


#----------------------------------------------------------------------

#==
```
_symm_tens_decomp(X, l, F0, rescaling, use_kernel, ::Val{:real})
```
Perform symmetric tensor decomposition with real weights, formulating the signed measure
as a difference of two positive measures ``\\mu - \\mu_-``.

==#
function _symm_tens_decomp(X, l, F0, rescaling, use_kernel, ::Val{:real})
    @info("Performing decomposition with real weights...")
    Y = vec(X[2:end])
    d = maxdegree(F0)
    s0 = _series_from_hpol(F0, X[1], Y, rescaling)

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

    ## Scale back the y and z variables ##
    X1[2:length(X),:] = X1[2:length(X),:].*rescaling
    X2[2:length(X),:] = X2[2:length(X),:].*rescaling


    ## Concatenate positive and negative components ##
    Xi = hcat(X1, X2)
    w = vcat(w1, -w2)

    return w, Xi
end


#----------------------------------------------------------------------
#==
```
_symm_tens_decomp(X, l, F0, rescaling, use_kernel, ::Val{:positive})
```
Perform symmetric tensor decomposition with strictly positive weights by formulating the problem
as a single positive measure ``\\mu``.
==#
function _symm_tens_decomp(X, l, F0, rescaling, use_kernel, ::Val{:positive})
    @info("Performing decomposition with positive weights...")
    Y = vec(X[2:end])
    d = maxdegree(F0)
    s0 = _series_from_hpol(F0, X[1], Y, rescaling)

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
   
    ## Scale back the y and z variables ##
    Xi[2:length(X),:] = Xi[2:length(X),:].*rescaling

    return w, Xi
end

# --- Public function ---

#----------------------------------------------------------------------
"""
```
w, Xi = symm_tens_decomp(X, l, F; rescaling=1, use_kernel=true, weight_type=:real)
```
Decompose a symmetric tensor represented as a homogeneous polynomial `F` into a weighted sum of rank-1 components.
Given a homogeneous polynomial `F0` of degree `d`, find weights `w` and points `Xi` such that
``F \\approx \\sum_{i} w_i \\langle X_i, x \\rangle^d``.

# Arguments
- `X`: Vector of polynomial variables
- `l`: Relaxation order for the moment problem
- `F0`: Homogeneous polynomial to decompose
- `rescaling`: Scaling factor to be set so decomposition points lie in radius 1 sphere (default: 1)
- `use_kernel`: Include kernel constraints from Hankel nullspace (default: true)
- `weight_type`: `:real` for signed weights or `:positive` for positive-only weights (default: `:real`)

# Returns
- `w`: Vector of weights
- `Xi`: Matrix of decomposition points (each column is a point)

"""
function symm_tens_decomp(X, l, F0; rescaling=1, use_kernel=true, weight_type=:real)
    return _symm_tens_decomp(X, l, F0, rescaling, use_kernel, Val(weight_type))
end

end #module


