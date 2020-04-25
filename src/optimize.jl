#using MathOptInterface

export getseries, getminimizers, getmeasure,
    minimize, maximize, optimize, minimize_ncl, minimize_tv

#----------------------------------------------------------------------
function getseries(M::MOM.Model)
    [series([M[:monomials][i]=>JuMP.value(M[:moments][i,k])
             for i in 1:length(M[:monomials])])
     for k in 1:M[:nu]]
end

#----------------------------------------------------------------------
function getminimizers(M::MOM.Model)
    s = getseries(M)[1]
    w, Xi = decompose(s);
    Xi
end

#----------------------------------------------------------------------
function getmeasure(M::MOM.Model)
    s = getseries(M)    

    w, Pts = decompose(s[1]);
    for k in 2:M[:nu]
        c, Xi = decompose(s[k])
        w = vcat(w, c*(-1)^(k-1))
        Pts= hcat(Pts,Xi)
    end
    
    return w, Pts
end

#----------------------------------------------------------------------
function JuMP.optimize!(M::MOM.Model)
    JuMP.optimize!(M.model)
    if JuMP.has_values(M.model)
        return JuMP.objective_value(M.model), M
    else
        println("Solver status: ", JuMP.termination_status(M.model))
        return nothing, M
    end
end

function JuMP.objective_value(M::MOM.Model)
    JuMP.objective_value(M.model)
end

function JuMP.value(M::MOM.Model)
    JuMP.value.(M[:moments])
end

function JuMP.set_optimizer(M::MOM.Model, optimizer)
    JuMP.set_optimizer(M.model, optimizer)
end

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
function minimize(fct, Eq, Pos,  X, d::Int64, optimizer; kwargs...)
    M = MOM.Model(X, d, optimizer; kwargs...)
    constraint_unitmass(M)
    constraint_zero(M,Eq...)
    constraint_nneg(M,Pos...)
    objective(M, fct)
    return JuMP.optimize!(M)
end

#----------------------------------------------------------------------
function maximize(fct, Eq, Pos,  X, d::Int64, optimizer; kwargs...)
    M = MOM.Model(X, d, optimizer; kwargs...)
    constraint_unitmass(M)
    constraint_zero(M,Eq...)
    constraint_nneg(M,Pos...)
    objective(M, fct, "sup")
    return JuMP.optimize!(M)
end

#----------------------------------------------------------------------
function optimize(C::Vector, X, d::Int64, optimizer; kwargs...)
    M  = MOM.Model(X,d;kwargs...)
    constraint_unitmass(M)
    for c in C
        if c[2] == "inf"
            objective(M,c[1], "inf")
        elseif c[2] == "sup"
            objective(M,c[1], "sup")
        elseif c[2] == "=0"
            constraint_zero(M,c[1])
        elseif c[2] == ">=0"
            constraint_nneg(M,c[1])
        elseif c[2] == "<=0"
            constraint_nneg(M,-c[1])
        end
    end
    set_optimizer(M,optimizer)
    return JuMP.optimize!(M)
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
