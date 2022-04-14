export constraint_zero, constraint_nneg, constraint_unitmass, constraint_moments

#----------------------------------------------------------------------
function add_constraint_zero(M::JuMP.Model, eq)
    p = eq*one(Polynomial{true,Float64})
    X = M[:variables]
    L = monomials(X,seq(0:2*M[:degree]-maxdegree(p)))
    for mn in L
        q = p*mn*one(Float64)
        for k in 1:M[:nu]
            @constraint(M,
                        sum(t.α*M[:moments][k,M[:index][t.x]] for t in q) == 0)
        end
    end
end

function add_constraint_zero(M::JuMP.Model, idx::Vector{Int64}, Veq)
    P = [ e*one(Polynomial{true,Float64}) for eq in Veq ]
    X = M[:variables]
    d0 = max([maxdegree(e) for e in P]...)
    L = monomials(X,seq(0:2*M[:degree]-d0))
    for mn in L
        @constraint(M,
                    sum(sum(t.α*M[:moments][idx[k],M[:index][t.x]] for t in P[k]*mn) for k in 1:length(idx)) == 0)
    end
end

#----------------------------------------------------------------------
function add_constraint_nneg(M::JuMP.Model, pos)
    p = pos*one(Polynomial{true,Float64})
    d0 = Int(ceil(maxdegree(p)/2))
    X = M[:variables]
    L = monomials(X, seq(0:M[:degree] - d0))
    N = length(L)
    if N == 1
        for k in 1:M[:nu]
            @constraint(M,
                        sum(t.α*M[:moments][k,M[:index][t.x]] for t in p)>=0)
        end
    else
        for k in 1:M[:nu]
            P = [ sum(t.α*M[:moments][k,M[:index][t.x*L[i]*L[j]]]
                      for t in p)+0
                  for i in 1:N, j in 1:N ]
            @constraint(M, Symmetric(P) in PSDCone())
        end
    end
end

function add_constraint_nneg(M::JuMP.Model, idx::Vector{Int64}, Veq)
    p = [ eq*one(Polynomial{true,Float64}) for eq in Veq ]
    X = M[:variables]
    d0 = Int(ceil(max([maxdegree(e) for e in p]...)/2))
    L = monomials(X, seq(0:M[:degree] - d0))
    N = length(L)
    if N == 1
        @constraint(M,
                    sum(sum(t.α*M[:moments][idx[k],M[:index][t.x]] for t in p[k]) for k in 1:length(idx)) >=0)
    else
        P = [ sum(sum(t.α*M[:moments][idx[k],M[:index][t.x*L[i]*L[j]]]
                      for t in p[k]) for k in 1:length(idx))+0
              for i in 1:N, j in 1:N ]
        @constraint(M, Symmetric(P) in PSDCone())
    end
end

#----------------------------------------------------------------------
function add_constraint_moment(M::JuMP.Model, v , p)
    for k in 1:M[:nu]
        @constraint(M,
                    sum(sum(t.α*M[:moments][k,M[:index][t.x]] for t in p) for k in 1:M[:nu])-v == 0)
    end
end

function add_constraint_moment(M::JuMP.Model, v, p::Vector)
    @constraint(M,
                sum(sum(t.α*M[:moments][k,M[:index][t.x]] for t in p[k]) for k in 1:length(p)) - v ==0)
end

function add_constraint_moment(M::JuMP.Model, v, idx::Vector{Int64}, p::Vector)
    @constraint(M,
                sum(sum(t.α*M[:moments][idx[k],M[:index][t.x]] for t in p[k]) for k in 1:length(idx)) - v ==0)
end


#----------------------------------------------------------------------
"""
```
constraint_zero(M, eqs ...)
```
Add to the moment program `M`, the constraints ``e \\star \\mu_i = 0`` for i in `1:M[:nu]` and e in eqs.

"""
function constraint_zero(M::JuMP.Model, eqs...)
    for e in eqs
        MOM.add_constraint_zero(M,e)
    end
end

"""
```
constraint_zero(M, i::Int, eqs ...)
```
Add to the moment program `M`, the constraint ``e \\star \\mu_i = 0`` for e in eqs.

"""
function constraint_zero(M::JuMP.Model, idx::Int64, eqs...)
    for e in eqs
        MOM.add_constraint_zero(M,[idx],[e])
    end
end


"""
```
constraint_zero(M, Eqs ...)
```
Add to the moment program `M`, the constraints ``\\sum_{i=1}^{\\nu} \\mathbf{e}_{i} \\star \\mu_i =0`` 
for ``\\mathbf{e}`` in Eqs where Eqs is a sequence of vectors of ``\\nu``=`M[:nu]` polynomials.

"""
function constraint_zero(M::JuMP.Model, eqs::Vector...)
    for e in eqs
        MOM.add_constraint_zero(M, collect(1:M[:nu]), e)
    end
end

function constraint_zero(M::JuMP.Model, idx::Vector{Int64}, eqs::Vector...)
    for e in eqs
        MOM.add_constraint_zero(M,idx,e)
    end
end

#----------------------------------------------------------------------
"""
```
constraint_nneg(M, eqs ...)
```
Add to the moment program `M`, the constraints ``e \\star \\mu_i \\succeq 0`` for i in `1:M[:nu]` and e in eqs.

"""
function constraint_nneg(M::JuMP.Model, eqs...)
    for e in eqs
        MOM.add_constraint_nneg(M,e)
    end
end

"""
```
constraint_nneg(M, i::Int64, eqs ...)
```
Add to the moment program `M`, the constraints ``e \\star \\mu_i \\succeq 0`` for  e in eqs.

"""
function constraint_nneg(M::JuMP.Model, idx::Int64, eqs...)
    for e in eqs
        MOM.add_constraint_nneg(M,[idx],[e])
    end
end

"""
```
constraint_nneg(M, Eqs ...)
```
Add to the moment program `M`, the constraints ``\\sum_{i=1}^{\\nu} \\mathbf{e}_{i} \\star \\mu_i \\succeq 0`` 
for ``\\mathbf{e}`` in Eqs where Eqs is a sequence of vectors of
``\\nu``=`M[:nu]` polynomials.

"""
function constraint_nneg(M::JuMP.Model, eqs::Vector...)
    for e in eqs
        MOM.add_constraint_nneg(M, collect(1:M[:nu]), e)
    end
end

function constraint_nneg(M::JuMP.Model, idx::Vector{Int64}, eqs::Vector...)
    for e in eqs
        MOM.add_constraint_nneg(M,idx,e)
    end
end

#----------------------------------------------------------------------
"""
```
constraint_moments(M, [m => c, ...])
```
Add to the moment program `M`, the constraints ``\\langle \\mu_i, m \\rangle - c = 0`` 
for i in `1:M[:nu]` and all pairs m=>c.

"""
function constraint_moments(M::JuMP.Model, moments::Vector)
    for c in moments
        MOM.add_constraint_moment(M, c[2], c[1]*one(Polynomial{true,Float64}))
    end
end

"""
```
constraint_moments(M, [m => c, ...], p)
```
Add to the moment program `M`, the constraints ``\\sum_{i=1}^{\\nu}\\langle {p}_{i} \\star \\mu_i, m\\rangle - c = 0`` 
for all pairs m=>c, where `p` is a vector of ``\\nu``=`M[:nu]` polynomials.

"""
function constraint_moments(M::JuMP.Model, moments::Vector, p::Vector)
    for c in moments
        MOM.add_constraint_moment(M, c[2], p*c[1]*one(Polynomial{true,Float64}))
    end
end

function constraint_moments(M::JuMP.Model, moments::Vector, idx::Vector{Int64}, P::Vector)
    for c in moments
        MOM.add_constraint_moment(M, c[2], idx, P*c[1]*one(Polynomial{true,Float64}))
    end
end

#----------------------------------------------------------------------
"""
```
constraint_unitmass(M)
```
Add to the moment program `M`, the constraints ``\\langle \\mu_i, 1 \\rangle - 1 = 0``  for i in `1:M[:nu]`.

"""

function constraint_unitmass(M::JuMP.Model)
    MOM.add_constraint_moment(M, 1, one(Polynomial{true,Float64}))
end

"""
```
constraint_unitmass(M, i::Int64, p)
```
Add to the moment program `M`, the constraints ``\\langle {p} \\star \\mu_i, 1 \\rangle - 1 = 0`` where `p` is a polynomial.

"""
function constraint_unitmass(M::JuMP.Model, idx::Int64, p)
    MOM.add_constraint_moment(M, 1, [idx], [p])
end

"""
```
constraint_unitmass(M, p)
```
Add to the moment program `M`, the constraints ``\\langle {p}_i \\star \\mu_i, 1 \\rangle - 1 = 0`` where `p` is a vector of ``\\nu``=`M[:nu]` polynomials.

"""
function constraint_unitmass(M::JuMP.Model, P::Vector)
    MOM.add_constraint_moment(M, 1, idx, P)
end

function constraint_unitmass(M::JuMP.Model, idx::Vector{Int64}, P::Vector)
    MOM.add_constraint_moment(M, 1, idx, P)
end
#----------------------------------------------------------------------
