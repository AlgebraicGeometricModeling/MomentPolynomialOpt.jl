export minimize

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
    m = MOM.Model(fct, Eq, Pos, X, d, optimizer)
    JuMP.optimize!(m.model)
     if JuMP.has_values(m.model)
        return JuMP.objective_value(m.model), m
    else
        println("Solver status: ", JuMP.termination_status(m.model))
        return nothing, m
    end
end

function minimize(pop,  X, d::Int64, optimizer)
    m = MOM.Model(pop,X, d, optimizer)
    JuMP.optimize!(m.model)
    if JuMP.has_values(m.model)
        return JuMP.objective_value(m.model), m
    else
        println("Solver status: ", JuMP.termination_status(m.model))
        return nothing, m
    end
end


#----------------------------------------------------------------------
# Minimize the nuclear norm, using a SDP matrix of size 2*N
function minimize_ncl(X, d, sigma, optimizer)

    m = MUS.Model(X, d, optimizer)
    JuMP.@objective(m.model, Min, sum( P[i,i] for i in 1:2*N) )

    if JuMP.has_values(m.model)
        return JuMP.objective_value(m.model), m
    else
        println("Solver status: ", JuMP.termination_status(m.model))
        return nothing, m
    end
end


function minimize_ncl(X, d, sigma)

    L = monoms(X,d)
    N = length(L)
    n = length(variables(sigma))

    mdl = Model(solver = CSDPSolver())

    @variable(mdl, P[1:2*N,1:2*N], SDP)

    I = Dict{DynamicPolynomials.Monomial{true}, Tuple{Int64,Int64}}()
    for i in 1:N
        for j  in 1:N
            m = L[i]*L[j]
	        if !haskey(I, m)
	            I[m] = (i,j)
	        else
                id = I[m]
	            if (i != id[2]) || (j != id[1])
                    @constraint(mdl, P[i,N+j]-P[id[1], N+id[2]] == 0)
	            end
	        end
	        if  degree(m) <= maxdegree(sigma) #cf != 0
                cf = sigma[m]
                @constraint(mdl, P[i,N+j]-cf == 0)
	        end
        end
    end

    #println(mdl)

    @objective(mdl, Min, sum( P[i,i] for i in 1:2*N) )

    status = solve(mdl)
    println("Objective value: ", getobjectivevalue(mdl)/2)
    R = getvalue(P)
    R[1:N,N+1:2*N]
end
