export minimize, minimize_ncl, minimize_tv

#------------------------------------------------------------------------
"""
```julia
v, m = minimize(f, L, [e1, e2, ...], [p1, p2, ...])
```
Compute the infimum of the functional associated to the polynomial f on the moments of the monomials in L*L such that the Hankel matrix with row and columns indexed by L is Semi_definite Positive and satisfying the equality constraints ei == 0 and the sign constraints pi >= 0.

f, ei, p1 should be polynomials in the variables X.

If the problem is feasible and has minimizers, it outputs
  - v: the minimum value 
  - m: the moment model of type MOM.Model 

Example
-------
```julia
using MomentTools

X  = @polyvar x1 x2
e1 = x1^2-2
e2 = (x2^2-3)*(x1*x2-2)
p1 = x1
p2 = 2-x2
v, m = minimize(-x1, [e1, e2], [p1, p2], X, 3)
```
The minimizers can be recovered as follows:
```julia
getminimizer(m)

[1.41421 1.73205; 1.41421 1.41421; 1.41421 -1.73205]
```
"""
function minimize(fct, Eq, Pos,  X, d::Int64, optimizer)
    M = MOM.Model(fct, Eq, Pos, X, d, optimizer)
    JuMP.optimize!(M)
    if JuMP.has_values(M.model)
        return JuMP.objective_value(M.model), M
    else
        println("Solver status: ", JuMP.termination_status(M.model))
        return nothing, M
    end
end

function minimize(pop,  X, d::Int64, optimizer)
    M = MOM.Model(pop,X, d, optimizer)
    JuMP.optimize!(M.model)
    if JuMP.has_values(M.model)
        return JuMP.objective_value(M.model), M
    else
        println("Solver status: ", JuMP.termination_status(M.model))
        return nothing, M
    end
end

#----------------------------------------------------------------------
# Minimize the nuclear norm, using a SDP matrix of size 2*N
function minimize_ncl(X, d, sigma, optimizer)

    M = MUS.Model(X, d, optimizer)
    JuMP.@objective(M.model, Min, sum( P[i,i] for i in 1:2*N) )

    if JuMP.has_values(M.model)
        return JuMP.objective_value(M.model), M
    else
        println("Solver status: ", JuMP.termination_status(M.model))
        return nothing, M
    end
end


#----------------------------------------------------------------------
# Minimize the nuclear norm, using a SDP matrix of size 2*N
function minimize_tv(X, d, sigma, optimizer)

    M = MUS.Model(X, d, optimizer)
    JuMP.@objective(M.model, Min, sum( P[i,i] for i in 1:2*N) )

    if JuMP.has_values(M.model)
        return JuMP.objective_value(M.model), M
    else
        println("Solver status: ", JuMP.termination_status(M.model))
        return nothing, M
    end
end

#----------------------------------------------------------------------
