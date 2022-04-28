export sos_decompose

using LinearAlgebra, JuMP, DynamicPolynomials

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



