export set_objective, set_objective_ncl, set_objective_tv

using JuMP

using DynamicPolynomials: coefficient

#----------------------------------------------------------------------
"""
```
set_objective(M, sense, p, mu)
```
Set the "inf" or "sup" objective function to  ``\\sum_{i=1}^{\\nu} \\langle p\\star \\mu_i, 1 \\rangle`` where `p` is a polynomial.

"""
function set_objective(M::JuMP.Model, sense::String, p, mu::Moments)
    if p != nothing
        obj = dot(mu, p)
        if sense == "inf"
            @objective(M, Min, obj)
        else
            @objective(M, Max, obj)
        end
    end
end


#----------------------------------------------------------------------
"""
```
set_objective_ncl(M)
```
Set the objective function of moment program to the nuclear norm or equivalently the trace of the moment matrices.
"""
function set_objective_ncl(M::JuMP.Model, mu)
    d = maxdegree(mu.basis)
    B = mu.basis[findall(x-> maxdegree(x) < div(d,2), mu.basis)]
    @objective(M, Min,  sum( dot(mu, b^2) for b in B))
end
#----------------------------------------------------------------------
"""
```
set_objective_tv(M)
```
Set the objective function of the moment program to the total variation of the moment sequence ``\\mu``, that is 
the sum  of the unit mass of the positive moment sequences ``\\sum_{i=1}^{\\nu} \\langle \\mu_i, 1\\rangle``.
"""
function set_objective_tv(M::JuMP.Model, mu)
    @objective(M, Min, dot(mu,1))
end

