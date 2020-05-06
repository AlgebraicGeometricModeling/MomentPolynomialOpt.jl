export constraint_zero, constraint_nneg, constraint_unitmass, constraint_moments

#----------------------------------------------------------------------
"""
```
constraint_zero(M, eqs ...)
```
Add to the moment program `M`, the constraints ``e \\star \\mu_i = 0`` for i in `1:M[:nu]` and e in eqs.

"""
function constraint_zero(M::MOM.Model, eqs...)
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
function constraint_zero(M::MOM.Model, idx::Int64, eqs...)
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
function constraint_zero(M::MOM.Model, eqs::Vector...)
    for e in eqs
        MOM.add_constraint_zero(M, collect(1:M[:nu]), e)
    end
end

function constraint_zero(M::MOM.Model, idx::Vector{Int64}, eqs::Vector...)
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
function constraint_nneg(M::MOM.Model, eqs...)
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
function constraint_nneg(M::MOM.Model, idx::Int64, eqs...)
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
function constraint_nneg(M::MOM.Model, eqs::Vector...)
    for e in eqs
        MOM.add_constraint_nneg(M, collect(1:M[:nu]), e)
    end
end

function constraint_nneg(M::MOM.Model, idx::Vector{Int64}, eqs::Vector...)
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
function constraint_moments(M::MOM.Model, moments::Vector)
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
function constraint_moments(M::MOM.Model, moments::Vector, p::Vector)
    for c in moments
        MOM.add_constraint_moment(M, c[2], p*c[1]*one(Polynomial{true,Float64}))
    end
end

function constraint_moments(M::MOM.Model, moments::Vector, idx::Vector{Int64}, P::Vector)
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

function constraint_unitmass(M::MOM.Model)
    MOM.add_constraint_moment(M, 1, one(Polynomial{true,Float64}))
end

"""
```
constraint_unitmass(M, i::Int64, p)
```
Add to the moment program `M`, the constraints ``\\langle {p} \\star \\mu_i, 1 \\rangle - 1 = 0`` where `p` is a polynomial.

"""
function constraint_unitmass(M::MOM.Model, idx::Int64, p)
    MOM.add_constraint_moment(M, 1, [idx], [p])
end

"""
```
constraint_unitmass(M, p)
```
Add to the moment program `M`, the constraints ``\\langle {p}_i \\star \\mu_i, 1 \\rangle - 1 = 0`` where `p` is a vector of ``\\nu``=`M[:nu]` polynomials.

"""
function constraint_unitmass(M::MOM.Model, P::Vector)
    MOM.add_constraint_moment(M, 1, idx, P)
end

function constraint_unitmass(M::MOM.Model, idx::Vector{Int64}, P::Vector)
    MOM.add_constraint_moment(M, 1, idx, P)
end
#----------------------------------------------------------------------
