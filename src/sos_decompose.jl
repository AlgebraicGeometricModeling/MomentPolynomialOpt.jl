export sos_decompose, rational_round,  w_cholesky, wsos, sym_matrix

using LinearAlgebra, JuMP, DynamicPolynomials


"""
       w, L = w_cholesky(A)

Compute the weighted Cholesky factorisation so that 'L^t*diagm(w)*L = A`

"""
function w_cholesky(A)
    s = size(A,1)
    L = fill(zero(A[1,1]),s,s)
    w = fill(zero(A[1,1]),s)
    for j in 1:s
        L[j,j] = 1
        w[j] = A[j,j] - sum(L[j,k]^2*w[k] for k in 1:j-1; init = zero(A[1,1]))
        for i in j+1:s
            L[i,j] = (A[i,j] - sum(L[i,k]*L[j,k]*w[k] for k in 1:j-1; init = zero(A[1,1])))/w[j]
        end
    end
    return w, L
end
#----------------------------------------------------------------------


function sym_matrix(pol,L)
    s = length(L)
    C =  coefficients(pol)
    M =  monomials(pol)
    S  = fill(zero(C[1]),s,s)
    Sc = Dict{typeof(M[1]), Int}()
    for i in 1:s
        for j in 1:s
            m = L[i]*L[j]
            if get(Sc,m,0) != 0
                Sc[m]+= 1
            else
                Sc[m] = 1
            end
            k = findfirst(x-> x == m, M)
            if k != nothing
                S[i,j] = C[k] 
            end
        end
    end
    
    for i in 1:s
        for j in 1:s
            m = L[i]*L[j]
            S[i,j] /= Sc[m]
        end
    end
    
    return S
end
#----------------------------------------------------------------------
function wsos(WS)
    sum(WS[1][i]*(WS[2][i])^2 for i in 1:length(WS[1]))
end

#----------------------------------------------------------------------
function rational_round(x::Number, k)
    BigInt(round(x*10^k; digits = 0))//10^k
end


function rational_round(p::Polynomial, k)
    dot(rational_round.(p.a,k),p.x)
end


"""
Decompose f in the truncated quadratic module associated to the constraints H=0, G>=0:

```
 WS, P, v, M = sos_decompose(f, H, G, X, d, optimizer; exact = false, round = Inf64 )
```
such that 

``  f = \\sum_k \\omega_{0,k} q_{0,k}^2 + \\sum_j (\\sum_k \\omega_{j,k} q_{j,k}^2)*G[j] + \\sum_i P[i]*H[i]              (+)  ``

where

  - f is the polynomial to decompose
  - H = [...] are the equality constraints
  - G = [...] are the non-negativity constraints
  - X is the set of variables
  - d is the order of the relaxation
 - `optimizer` (optional) is the optimizer used to solve the SDP program. The default value is `MPO[:optimizer]`
 - if `exact = true` then an exact decomposition with rational coefficients is computed.
 - if the option `round = p` is provided, then a round of at most `p`digits is used.

In the output decomposition, 

 - WS[1] = [w0, Q0] where w0 is an array of weights ``\\omega_{0,k}`` in (+) and Q0 is an array of polynomials (q_{0,k} in (+)) of degree ``\\le 2d -\\deg(G[i])``
 - WS[j+1] as WS[1] for the weights ``\\omega_{j,k}`` and polynomials ``q_j,k``in (+)
 - P[i] is a polynomial of degree ``\\le 2d - \\deg(H[i])``
 - v is the maximal smallest eigenvalue of S0 such that ``(X^d)^t S0 (X^d) = \\sum_k \\omega_{0,k} q_{0,k}^2`` in (*)
 - `M`is the JuMP optimization model 

"""
function sos_decompose(f,  H::AbstractVector, G::AbstractVector, X, d::Int64, optimizer = MPO[:optimizer]; exact = false, round = Inf64)


    M = MaxEigenModel(f, H, G, X, d, optimizer)
    
    optimize!(M)
    lambda = objective_value(M)

    if lambda <0
        @warn "not in the interior of the truncated quadratic module"
        println(">  ", eigen(value.(M[:Q0])).values)
        #        return nothing
    end
    
    L0 = M[:monomials]
    MP = M[:mP]
    MQ = M[:mQ]
    
    #Q0 = L0'*value.(M[:Q0])*L0



    if length(H) >0 
        P = [ dot(MP[i],value.(M[:P][i])) for i in 1:length(H)]
    else
        P = []
    end

    #for i in 1:length(G) println(">>  ", eigen(value.(M[:Q][i])).values) end

    if exact != true
    #== Approximate decomposition ==#
        WS = fill([], length(G)+1)
        w0, Lw0 = w_cholesky(value.(M[:Q0]))
        WS[1] = [w0, Lw0'*L0]
        for i in 1:length(G)
            w, Lw = w_cholesky(value.(M[:Q][i]))
            WS[i+1] = [w, Lw' * MQ[i]]
        end

        return WS, P, lambda, M
    end

    #== Exact decomposition ==#
        
    k = max(Int(ceil(log10(1/lambda)))+3, 0)
    if round < k
        k = Int(round)
    end

    @info "rounding $k digit(s)"
    
    WS = fill([], length(G)+1)
    for i in 1:length(G)
        w, Lw = w_cholesky(value.(M[:Q][i]))
        Lk = rational_round.(Lw,k)
        wk = rational_round.(w,k)
        println(" >>>> ", norm(Lw-Float64.(Lk)))
        println(" >>>> ", norm(value.(M[:Q][i]) - Float64.(Lk*diagm(wk)*Lk')))
        WS[i+1] = [wk, Lk' * MQ[i]]
    end


    q0e = dot(Rational{BigInt}.(coefficients(f)), monomials(f))

    if length(H)>0
        Pe = rational_round.(P,k)
        q0e = q0e  - dot(H,Pe) 
    else
        Pe = []
    end

    if length(G)>0
        q0e = q0e - dot(G, wsos.(WS[2:end]))
    end

    Q0r = rational_round.(value.(M[:Q0]),k)

    q0e = q0e - L0'*Q0r*L0
    
    QS0 = Q0r + sym_matrix(q0e, L0)

    w0, Lw0 = w_cholesky(QS0)

    WS[1] = [w0, Lw0'*L0]
    
    return WS, Pe, lambda, M
end

function sos_decompose(f, X = variables(f), d::Int64 = Int64(ceil(maxdegree(f)/2)),
                       optimizer = MPO[:optimizer]; exact = false, round = Inf64)

    WS0, WS, P, v, M = sos_decompose(f, [], [], X, d; exact, round)
    return WS0, v, M

end

#----------------------------------------------------------------------
function m_sqrt(m::DynamicPolynomials.Monomial)
    return prod(m.vars.^div.(m.z,2))
end

function m_div(m::DynamicPolynomials.Monomial, m1::DynamicPolynomials.Monomial)
    return prod(m.vars.^(m.z-m1.z))
end

function is_null(f::DynamicPolynomials.Polynomial)
    return (length(f.x)==0)
end

#----------------------------------------------------------------------

