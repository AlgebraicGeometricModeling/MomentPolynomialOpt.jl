export constraint_zero, constraint_nneg, constraint_unitmass, constraint_moments

using LinearAlgebra

using DynamicPolynomials: coefficient

#----------------------------------------------------------------------
"""
```
add_constraint_zero(M, eq, mu::Moments)
```
Add to the moment program `M`, the constraints ``eq\\star \\mu = 0``, that is all the constraints ``\\mu(p*m)==0`` for ``degree(m) \\le 2*degree(M)-degree(p)``.

"""
function add_constraint_zero(M::JuMP.Model,  mu::Moments, eq)
    p = convert_Float64(eq)
    X = variables(mu.basis) 
    L = monomials(X,0:maxdegree(mu.basis) - maxdegree(p))
    for mn in L
        @constraint(M, dot(mu, p*mn) == 0)
    end
end

#----------------------------------------------------------------------
"""
```
add_constraint_nneg(M, mu::Moments, g)
```
Add to the moment program `M`, the constraints ``g \\star \\mu \\succeq 0``, that is the moment matrix of ``g\\star \\mu`` is PSD.

"""
function add_constraint_nneg(M::JuMP.Model, mu::Moments, g)
    p = convert_Float64(g)
    d0 = Int(ceil(maxdegree(p)/2))
    d = div(maxdegree(mu.basis),2)
    X = variables(mu.basis) #X = M[:variables]
    L = monomials(X, 0:d-d0)
    N = length(L)
    if N == 1
        @constraint(M, dot(mu, p) >= 0)
    else
        P = [ dot(mu, p*L[i]*L[j]) for i in 1:N, j in 1:N ]
        @constraint(M, Symmetric(P) in PSDCone())
    end
end

function add_constraint_nneg(M::JuMP.Model, mu::Moments)
    p = convert_Float64(mu.basis[1])
    d = Int(floor(maxdegree(mu.basis)/2))
    X = variables(mu.basis)
    L = monomials(X, 0:d)
    N = length(L)
    P = [ dot(mu, p*L[i]*L[j]) for i in 1:N, j in 1:N ]
    @constraint(M, Symmetric(P) in PSDCone())
end

#----------------------------------------------------------------------
"""
```
add_constraint_moment(M, mu::Moments,  p => v)
```
Add to the moment program `M`, the constraint ``\\mu(p) - v == 0``.

"""
function add_constraint_moment(M::JuMP.Model, mu:: Moments, pv ::Pair )
    @constraint(M,
                sum(coefficient(t)*dot(mu,monomial(t)) for t in terms(pv[1])) -pv[2] == 0)
end

#----------------------------------------------------------------------
"""
```
add_constraint_unitmass(M, mu)
```
Add to the moment program `M`, the constraint ``\\langle \\mu, 1 \\rangle - 1 = 0``.

"""
function add_constraint_unitmass(M::JuMP.Model, mu::Moments)
    @constraint(M, mu.values[1] - 1 == 0)
end
#----------------------------------------------------------------------
