export set_objective, set_objective_ncl, set_objective_tv

using JuMP

using DynamicPolynomials: coefficient

#----------------------------------------------------------------------
"""
```
set_objective(M, sense, p)
```
Set the "inf" or "sup" objective function to  ``\\sum_{i=1}^{\\nu} \\langle p\\star \\mu_i, 1 \\rangle`` where `p` is a polynomial.

"""
function set_objective(M::JuMP.Model, sense::String, p)
    if p != nothing
        f = convert_Float64(p)
        obj = sum(sum(coefficient(t)*M[:moments][k,M[:index][monomial(t)]] for t in f) for k in 1:M[:nu])
        if sense == "inf"
            @objective(M, Min, obj)
        else
            @objective(M, Max, obj)
        end
    end
end


"""
```
set_objective(M, sense, p)
```
Set the "inf" or "sup" objective function to  ``\\sum_{i=1}^{\\nu} \\langle p_{i} \\star \\mu_i, 1 \\rangle`` where `p` is a vector of ``\\nu``=`M[:nu]` polynomials. 

"""
function set_objective(M::JuMP.Model, sense, idx::Vector{Int64}, p::Vector)
    if p != nothing
        f = convert_Float64(p)
        obj = sum(sum(coefficient(t)*M[:moments][idx[k],M[:index][monomial(t)]] for t in f[k]) for k in 1:length(idx))
        if sense == "inf"
            @objective(M, Min, obj)
        else
            @objective(M, Max, obj)
        end
    end
end


function set_objective(M::JuMP.Model, sense, pol::Vector)
    if pol != nothing
        set_objective(M, collect(1:M[:nu]), pol, sense)
    end
end


"""
```
set_objective(M, sense, i, p)
```
Set the "inf" or "sup" objective function to  `` \\langle p \\star \\mu_i, 1 \\rangle`` where `p` is a polynomial. 
"""
function set_objective(M::JuMP.Model, sense::String, pol, idx::Int64)
    if pol != nothing
        set_objective(M, sense, [idx], [pol])
    end
end


# """
#  Add as objective function the linear functional associated to the polynomial pol to minimize.
# """
# function JuMP.set_objective(M::MOM.Model, sense, idx::Vector{Int64}, pol::Vector)
#     if pol != nothing
#         f = [p*one(Polynomial{true,Float64}) for p in pol]
#         MOM.set_objective(M, sense, idx, f)
#     end
# end
#----------------------------------------------------------------------
"""
```
set_objective_ncl(M)
```
Set the objective function of moment program to the nuclear norm or equivalently the sum of 
the traces of the moment matrices.
"""
function set_objective_ncl(M)
    B = M[:basis]
    @objective(M, Min,  sum(sum(M[:moments][k,M[:index][B[i]^2]] for i in 1:length(B)) for k in 1:M[:nu]))
end
#----------------------------------------------------------------------
"""
```
set_objective_ncl(M)
```
Set the objective function of moment program to the total variation of the moment sequences ``\\mu_i``, that is 
the sum  of the unit mass of the positive moment sequences ``\\sum_{i=1}^{\\nu} \\langle \\mu_i, 1\\rangle``.
"""
function set_objective_tv(M::JuMP.Model)
    set_objective(M,"inf", collect(1:M[:nu]), ones(M[:nu]))
end

