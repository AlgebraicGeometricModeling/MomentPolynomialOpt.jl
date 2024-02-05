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
WS0, WS, P, v, M = sos_decompose(f, G, H, X, d, optimizer; exact = false, round = Inf64 )
```
such that 

``  f = \\sum_k \\omega_{0,k} q_{0,k}^2 + \\sum_j (\\sum_k \\omega_{j,k} q_{j,k}^2)*G[j] + \\sum_i P[i]*H[i]              (+)  ``

where

  - f is the polynomial to decompose
  - G = [...] are the non-negativity constraints
  - H = [...] are the equality constraints
  - X is the set of variables
  - d is the order of the relaxation
 - `optimizer` (optional) is the optimizer used to solve the SDP program. The default value is `MMT[:optimizer]`
 - if `exact = true` then an exact decomposition with rational coefficients is computed.
 - if the option `round = p` is provided, then a round of at most `p`digits is used.

In the output decomposition, 

 - WS0 = [w0, Q0] where w0 is an array of weights ``\\omega_{0,k}`` in (+) and Q0 is an array of polynomials (q_{0,k} in (+)) of degree ``\\le 2d -\\deg(G[i])``
- WS = [WS1, WS2, ...] with WSi as WS0 for the weights ``\\omega_{j,k}`` and polynomials ``q_j,k``in (+)
 - P[i] is a polynomial of degree ``\\le 2d - \\deg(H[i])``
 - v is the maximal smallest eigenvalue of S0 such that ``(X^d)^t S0 (X^d) = \\sum_k \\omega_{0,k} q_{0,k}^2`` in (*)
 - `M`is the JuMP optimization model 

"""
function sos_decompose(f,  G::AbstractVector, H::AbstractVector, X, d::Int64, optimizer = MMT[:optimizer]; exact = false, round = Inf64)


    M = MaxEigenModel(f, G, H, X, d, optimizer)
    
    optimize!(M)
    lambda = objective_value(M)

    if lambda <0
        @warn "not in the interior of the truncated quadratic module"
        #        return nothing
    end
    
    L0 = M[:monomials]
    MP = M[:mP]
    MQ = M[:mQ]
    
    #Q0 = L0'*value.(M[:Q0])*L0

    #println(">  ", eigen(value.(M[:Q0])).values)

    if length(H) >0 
        P = [ dot(MP[i],value.(M[:P][i])) for i in 1:length(H)]
    else
        P = []
    end

    #for i in 1:length(G) println(">>  ", eigen(value.(M[:Q][i])).values) end

    if exact != true
        WS = fill([], length(G))
        for i in 1:length(G)
            w, Lw = w_cholesky(value.(M[:Q][i]))
            WS[i] = [w, Lw' * MQ[i]]
        end
        w0, Lw0 = w_cholesky(value.(M[:Q0]))
        return [w0, Lw0'*L0], WS, P, lambda, M
    end

    #== Exact decomposition ==#
        
    k = max(Int(ceil(log10(1/lambda)))+3, 0)
    if round < k
        k = Int(round)
    end

    @info "rounding $k digit(s)"
    
    WS = fill([], length(G))
    for i in 1:length(G)
        w, Lw = w_cholesky(value.(M[:Q][i]))
        Lk = rational_round.(Lw,k)
        wk = rational_round.(w,k)
        #println(" >>>> ", norm(Lw-Float64.(Lk)))
        #println(" >>>> ", norm(value.(M[:Q][i]) - Float64.(Lk*diagm(wk)*Lk')))
        WS[i] = [wk, Lk' * MQ[i]]
    end


    Q0e = dot(Rational{BigInt}.(coefficients(f)), monomials(f))

    if length(H)>0
        Pe = rational_round.(P,k)
        Q0e = Q0e  - dot(H,Pe) 
    else
        Pe = []
    end

    if length(G)>0
        Q0e = Q0e - dot(G, wsos.(WS))
    end

    QS0 = sym_matrix(Q0e, L0)

    w0, Lw0 = w_cholesky(QS0)

    WS0 = [w0, Lw0'*L0]
    
    return WS0, WS, Pe, lambda, M
end

function sos_decompose(f, X = variables(f), d::Int64 = Int64(ceil(maxdegree(f)/2)),
                       optimizer = MMT[:optimizer]; exact = false, round = Inf64)

    WS0, WS, P, v, M = sos_decompose(f, [], [], X, d; exact, round)
    return WS0, v, M

end

#----------------------------------------------------------------------
function m_sqrt(m::Monomial)
    return prod(m.vars.^div.(m.z,2))
end

function m_div(m::Monomial, m1::Monomial)
    return prod(m.vars.^(m.z-m1.z))
end

function is_null(f::Polynomial)
    return (length(f.x)==0)
end

#----------------------------------------------------------------------

function Q_matrix(g, n)
    Q  = fill(zero(DynamicPolynomials.coefficients(g)[1]),n,n)
    for t in g
        i = maxdegree(t)
        c = DynamicPolynomials.coefficient(t)
        if i <= n-1
            if i == 0
                Q[1,1]= c
            else
                Q[1,i+1] = c/2
            end
        else
            if i == 2n-2
                Q[n,n]= c
            else
                Q[i+2-n,n] = c/2
            end
        end
    end
    return Symmetric(Q)
end

function H_matrices(n, d)
    [
        [(i+j==k+2 ? 1//1 : 0//1 ) for i in 1:n, j in 1:n ]
        for k in 0:d
    ]
end

function H_proj(Q, H)
    sum(h*tr(h*Q)/tr(h*h) for h in H) 
end


#----------------------------------------------------------------------
"""
Decompose g>0 at the real roots of f as a rational Sum of Squares and a multiple of f with rational coefficients.

  - Q rational matrix which is positive definite if g>0 at the real roots of f
  - q rational factor 

such that g = x^t Q x + q*f.

The options are:

   - rounding = (Int64)k upper bound on the precision used in the rounding step
   - optimizer = optimizer tool used for solving the SDP problem
   - verbose = true or false to specify if information during the computation is printed or not

"""
function esos_decompose(g, f; rounding = 10000, optimizer = MMT[:optimizer], verbose = false)
    n = maxdegree(f)
    dg = maxdegree(g)
    dq = n-dg

    X  = variables(f)
    L  = monomials(X,0:2n-2)
    B  = monomials(X,0:n-1)
    Lq = monomials(X,0:dq)

    cf = MultivariateSeries.matrixof([f],L)
    cg = MultivariateSeries.matrixof([g],L)

    M = JuMP.Model(optimizer)

    # The variables are lambda, the symmetric matrix Q and the coefficients of q
    @variable(M, lambda)
    # @variable(M, Q[1:n, 1:n], PSD)
    @variable(M, Q[1:n, 1:n], Symmetric)
    @variable(M, q[1:dq+1])


    # Positivity constraint: Q - lambda*I >=0
    @constraint(M, Q - diagm(fill(lambda,n)) in PSDCone())

    # The constraints g = x^t Q x + q*f
    for k in 1:length(L)
        s = min(n,k)
        @constraint(M,
                    sum(Q[s+1-i,k+i-s] for i in 1:min(k,2n-k))
                    + sum(q[i]*cf[k+1-i] for i in 1:min(dq,k))
                    - cg[1,k] == 0)
    end    

    # Objectiv function = lambda to maximize the smallest eigenvalue of Q

    @objective(M, Max, lambda)

    JuMP.optimize!(M)

    lambda = JuMP.value.(lambda)

    Qopt = JuMP.value.(Q)
    qopt = JuMP.value.(q)
    gopt = B'*Qopt*B + (qopt'*Lq)*f

    H = H_matrices(n,2n-2)

    popt = - (qopt'*Lq)*f + g
    Qpopt = Q_matrix(popt,n)

    epol = g - (qopt'*Lq)*f - B'*Qopt*B

    rho = norm(epol)
    if verbose 
        println("Minimal eigenvalue: ", lambda)
        #Ev = eigen(Qopt).values
        #println("Eigenvalues: ", Ev)
    end

    
    # Rounding
    #delta=0.99*(min(Ev...)-rho)/(n+ sqrt(n)*(n-1)*norm(f))
    delta=0.99*(lambda-rho)/(n+ sqrt(n)*(n-1)*norm(f))
    #println("delta: ", delta)

    k = Int64(ceil(log10(1/delta)))
    k = min(k, rounding)

    if verbose
        println("Error: ", rho)
        println("Delta: ", delta)
        println("Rounding precision: $k digit(s)")
    end

    Qround = Int64.(round.(Qopt*10^k; digits = 0))//10^k
    qround = Int64.(round.(qopt*10^k; digits = 0))//10^k
 
    
    # Projection
    qexact = (qround'*Lq)
    p  = -qexact*f + g
    Qp = Q_matrix(p,n)
    Qexact = Qround + H_proj(Qp-Qround, H)
    
    
    if verbose 
        # Certificate
        delta = B'*Qexact*B+qexact*f-g
        # println("delta: ",delta)
        @assert delta == 0
    end
    
    Qexact, qexact
end



