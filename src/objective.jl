export objective, objective_ncl, objective_tv

"""
```
objective(M, p, sense="inf")
```
Set the "inf" or "sup" objective function to  ``\\sum_{i=1}^{\\nu} \\langle p\\star \\mu_i, 1 \\rangle`` where `p` is a polynomial.

"""
function objective(M::MOM.Model, pol, sense="inf")
    if pol != nothing
        MOM.set_objective(M,pol,sense)
    end
end

"""
```
objective(M, p, sense="inf")
```
Set the "inf" or "sup" objective function to  ``\\sum_{i=1}^{\\nu} \\langle p_{i} \\star \\mu_i, 1 \\rangle`` where `p` is a vector of ``\\nu``=`M[:nu]` polynomials. 

"""
function objective(M::MOM.Model, pol::Vector, sense="inf")
    if pol != nothing
        MOM.set_objective(M, collect(1:M[:nu]), pol, sense)
    end
end


"""
```
objective(M, i, p, sense="inf")
```
Set the "inf" or "sup" objective function to  `` \\langle p \\star \\mu_i, 1 \\rangle`` where `p` is a polynomial. 
"""
function objective(M::MOM.Model, idx:: Int64, pol, sense="inf")
    if pol != nothing
        MOM.set_objective(M, [idx], [pol], sense)
    end
end


"""
 Add as objective function the linear functional associated to the polynomial pol to minimize.
"""
function objective(M::MOM.Model, idx::Vector{Int64}, pol::Vector, sense="inf")
    if pol != nothing
        f = [p*one(Polynomial{true,Float64}) for p in pol]
        MOM.set_objective(M,idx,f, sense)
    end
end
#----------------------------------------------------------------------
"""
```
objective_ncl(M)
```
Set the objective function of moment program to the nuclear norm or equivalently the sum of 
the traces of the moment matrices.
"""
function objective_ncl(M)
    B = M[:basis]
    @objective(M.model, Min,  sum(sum(M[:moments][k,M[:index][B[i]^2]] for i in 1:length(B)) for k in 1:M[:nu]))
end
#----------------------------------------------------------------------
"""
```
objective_ncl(M)
```
Set the objective function of moment program to the total variation of the moment sequences ``\\mu_i``, that is 
the sum  of the unit mass of the positive moment sequences ``\\sum_{i=1}^{\\nu} \\langle \\mu_i, 1\\rangle``.
"""
function objective_tv(M)
    objective(M,collect(1:M[:nu]), ones(M[:nu]),"inf")
end

