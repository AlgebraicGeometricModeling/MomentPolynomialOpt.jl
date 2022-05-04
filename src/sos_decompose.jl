export sos_decompose

using LinearAlgebra, JuMP, DynamicPolynomials, MultivariateSeries

"""
Decompose f in the truncated quadratic module associated to the constraints H=0, G>=0:

```
s, p, q, v, M = sos_decompose(f,H,G,X,d, optimizer)
```
such that 

``  f = s + \\sum_i p[i]*H[i]+ \\sum_j q[j]*G[j]      (*)  ``

where 
 - p[i] is a polynomial of degree ``\\le 2d -\\deg(H[i])``
 - q[i] is a Sum of Squares of degree ``\\le 2d -\\deg(G[i])``
 - s is a Sum of Squares if v >= 0.
 - v is the maximal smallest eigenvalue of Q such that s = (X^d)^t Q (X^d) in (*)
 - `optimizer` (optional) is the optimizer used to solve the SDP program. The default value is `MMT[:optimizer]`.

"""
function sos_decompose(f, H::AbstractVector, G::AbstractVector, X, d::Int64, optimizer = MMT[:optimizer])

    M = JuMP.Model(with_optimizer(optimizer))
    
    M[:type] = :polynomial
    
    Mh = [monomials(X,0:2*d-maxdegree(h)) for h in H]

    Lg = [monomials(X,0:d-Int64(ceil(maxdegree(g)/2))) for g in G]
    L0 = monomials(X,0:d)

    @variable(M, lambda)

    a = [@variable(M, [1:length(Mh[i])], base_name="a$i") for i in 1:length(H)]

    Q = [@variable(M, [1:length(Lg[i]), 1:length(Lg[i])], PSD, base_name="Q$i") for i in 1:length(G)]

    n = length(L0)
    Q0 = @variable(M, [1:n, 1:n], Symmetric, base_name="Q0") 

    # Positivity constraint: Q - lambda*I >=0
    @constraint(M, Q0 - diagm(fill(lambda,n)) in PSDCone())
        
    M[:a] = a
    M[:Q] = Q
    M[:Q0]= Q0

    P = f - L0'*Q0*L0 -
            sum( H[i]*sum(a[i][j]*Mh[i][j] for j in 1:length(Mh[i])) for i in 1:length(H)) -
            sum( G[i]*(Lg[i]'*Q[i]*Lg[i]) for i in 1:length(G) )


    for c in coefficients(P)
        @constraint(M, c == 0)
    end

    @objective(M, Max, lambda)

    
    optimize!(M)
    v = objective_value(M)

    s0 = L0'*value.(Q0)*L0

    h = [ dot(Mh[i],value.(a[i])) for i in 1:length(H)]
    q = [ Lg[i]'*value.(Q[i])*Lg[i] for i in 1:length(G)]
    return s0, h, q, v, M

end



#----------------------------------------------------------------------
export exact_decompose


function Q_matrix(g, n)
    Q  = fill(zero(coefficients(g)[1]),n,n)
    for t in g
        i = maxdegree(t)
        c = coefficient(t)
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
function exact_decompose(g, f; rounding = 10000, optimizer = MMT[:optimizer], verbose = false)
    n = maxdegree(f)
    dg = maxdegree(g)
    dq = n-dg

    X  = variables(f)
    L  = monomials(X,0:2n-2)
    B  = monomials(X,0:n-1)
    Lq = monomials(X,0:dq)

    cf = MultivariateSeries.matrixof([f],L)
    cg = MultivariateSeries.matrixof([g],L)

    M = JuMP.Model(with_optimizer(optimizer))

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



