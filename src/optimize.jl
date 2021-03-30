export get_series, get_minimizers, get_measure,
    minimize, maximize, optimize, minimize_ncl, minimize_tv

#----------------------------------------------------------------------
"""
```
get_series(M)
```
Return the vector of ``\\nu``=`M[:nu]` series of optimal moments of the optimized moment program `M`.
"""
function get_series(M::MOM.Model)
    MOM.get_series(M)
    [series([M[:monomials][i]=>JuMP.value(M[:moments][k,i])
             for i in 1:length(M[:monomials])])
     for k in 1:M[:nu]]
end

#----------------------------------------------------------------------
"""
```
get_minimizers(M)
```
Return the minimizer points  of the optimized moment program `M`.

```julia
get_minimizer(M)

[1.41421 1.73205; 1.41421 1.41421; 1.41421 -1.73205]
```
"""
function get_minimizers(M::MOM.Model)
    s = get_series(M)[1]
    w, Xi = MultivariateSeries.decompose(s);
    Xi
end

#----------------------------------------------------------------------
"""
```
w, Xi = get_measure(M, lambda = [(-1)^(k-1) for k in 1:M[:nu]])
```
Return the approximation of the moment sequence ``\\sum_{i=1}^{\\nu} \\lambda_i \\mu_i`` as weighted sum of Dirac measures: ``\\sum_{k=1}^{r} \\omega_k \\delta_{\\xi_k}`` where 

- `w` is the vector of weights of the Dirac measures.
- `Xi` is matrix of ``n\\times r`` support points of the corresponding Dirac measures. The column `Xi[:,i]` is the support point ``\\xi_{i}`` of the ith Dirac measure and its weights is `w[i]`.

```julia
w, Xi = get_measure(M)

([0.1541368146508854, 0.5889741915171074, 0.256888993597116], [1.4142135624216647 1.414213562080608 1.4142135620270329; -1.732052464639053 1.4141771454788292 1.7319839273833693])
```
"""
function get_measure(M::MOM.Model, lambda::Vector = [(-1)^(k-1) for k in 1:M[:nu]])
    s = get_series(M)    
    w, Pts = MultivariateSeries.decompose(s[1]);
    for k in 2:M[:nu]
        c, Xi = MultivariateSeries.decompose(s[k])
        w = vcat(w, c*lambda[k])
        Pts= hcat(Pts,Xi)
    end
    
    return w, Pts
end

function get_measure(M::MOM.Model, e::Float64, lambda::Vector = [(-1)^(k-1) for k in 1:M[:nu]])
    s = get_series(M)    
    w, Pts = MultivariateSeries.decompose(s[1], MultivariateSeries.eps_rkf(e));
    for k in 2:M[:nu]
        c, Xi = MultivariateSeries.decompose(s[k],MultivariateSeries.eps_rkf(e))
        w = vcat(w, c*lambda[k])
        Pts= hcat(Pts,Xi)
    end
    
    return w, Pts
end


    
#----------------------------------------------------------------------
"""
```julia
v, M = optimize(M)
```
Run the optimizer on the moment program `M` and output the objective_value `v` and the moment program `M`. If the optimization program has no value, it returns `nothing` and `M`.
"""
function optimize(M::MOM.Model)
    JuMP.optimize!(M.model)
    if JuMP.has_values(M.model)
        return JuMP.objective_value(M.model), M
    else
        println("Solver status: ", JuMP.termination_status(M.model))
        return nothing, M
    end
end

#----------------------------------------------------------------------
function JuMP.objective_value(M::MOM.Model)
    JuMP.objective_value(M.model)
end

function JuMP.value(M::MOM.Model)
    JuMP.value.(M[:moments])
end

"""
```julia
v, M = JuMP.set_optimizer(M, optimizer)
```
Set the optimizer of the moment program `M` to the dual optimizer of `optimizer`.
"""
function JuMP.set_optimizer(M::MOM.Model, optimizer)
    if M[:dual]
        set_optimizer(M.model, Dualization.dual_optimizer(optimizer))
    else
        set_optimizer(M.model, optimizer)
    end
end

#------------------------------------------------------------------------
"""
```julia
v, M = minimize(f, [e1, e2, ...], [p1, p2, ...], X, d, optimizer)
```
Compute the infimum of `f` under the constraints ``e_i`` =0 and ``p_i \\geq 0`` using a relaxation of order `d` on the moments in the variable `X` and the optimizer `M`.

``f, e_i, p_i `` should be polynomials in the variables X.

If the problem is feasible and has minimizers, it outputs
  - v: the minimum value 
  - M: the moment model of type MOM.Model 

Example
-------
```julia
using MomentTools

X  = @polyvar x1 x2
e1 = x1^2-2
e2 = (x2^2-3)*(x1*x2-2)
p1 = x1
p2 = 2-x2
v, M = minimize(-x1, [e1, e2], [p1, p2], X, 3)
```
To recover the optimal values, see [`get_minimizers`](@ref), [`get_measure`](@ref), [`get_series`](@ref).

"""
function minimize(fct, Eq::Vector, Pos::Vector,  X, d::Int64, optimizer; kwargs...)
    M = MOM.Model(X, d, optimizer; kwargs...)
    constraint_unitmass(M)
    constraint_zero(M,Eq...)
    constraint_nneg(M,Pos...)
    if fct != nothing
        objective(M, fct)
    else
        objective(M, one(Polynomial{true,Float64}))
    end
    JuMP.optimize!(M.model)
    v = objective_value(M.model)
    return v, M
end

#----------------------------------------------------------------------
"""
```julia
v, M = maximize(f, [e1, e2, ...], [p1, p2, ...], X, d, optimizer)
```
Similar to the function `minimize` but compute the supremun of `f`.
"""
function maximize(fct, Eq::Vector, Pos::Vector,  X, d::Int64, optimizer; kwargs...)
    M = MOM.Model(X, d, optimizer; kwargs...)
    constraint_unitmass(M)
    constraint_zero(M,Eq...)
    constraint_nneg(M,Pos...)
    if fct != nothing
        objective(M, 1, fct, "sup")
    else
        objective(M, 1, 1, "sup")
    end
    JuMP.optimize!(M.model)
    v = objective_value(M.model)
    return v, M
end

#----------------------------------------------------------------------
"""
```julia
v, M = optimize([(f, set), ...], X, d, optimizer)
```
Solve the moment program of relaxation of order `d` in the variables `X`, defined by the constraint or objective paires `(f, set)` where 
`f` is a polynomial and `set` is a string
   
- "inf", "sup" to define the objective.
- "=0" to define zero constraints:
- ">=0", "<=0" to define sign constraints

It outputs
  - v: the optimal value 
  - M: the optimized moment program of type MOM.Model 

Example
-------
```julia
using MomentTools

X  = @polyvar x1 x2
e1 = x1^2-2
e2 = (x2^2-3)*(x1*x2-2)
p1 = x1
p2 = 2-x2
v, M = optimize([(-x1, "inf"), (e1, "=0"), (e2, "=0"), (p1, ">=0"), (p2>=0)], X, 3)
```    

To recover the optimal values, see [`get_minimizers`](@ref), [`get_measure`](@ref), [`get_series`](@ref).

"""
function optimize(C::Vector, X, d::Int64, optimizer; kwargs...)
    M  = MOM.Model(X,d;kwargs...)
    constraint_unitmass(M)
    wobj = false
    for c in C
        if c[2] == "inf" || c[2] == "min"
            objective(M,c[1], "inf")
            wobj = true
        elseif c[2] == "sup" || c[2] == "max"
            objective(M,c[1], "sup")
            wobj = true
        elseif c[2] == "=0"
            constraint_zero(M,c[1])
        elseif c[2] == ">=0"
            constraint_nneg(M,c[1])
        elseif c[2] == "<=0"
            constraint_nneg(M,-c[1])
        elseif isa(c[2],AbstractVector)
            constraint_nneg(M,c[1]-c[2][1])
            constraint_nneg(M,-c[1]+c[2][2])
        end
    end
    if !wobj
        objective(M,one(C[1][1]), "sup")
    end
    set_optimizer(M,optimizer)
    JuMP.optimize!(M.model)
    v = objective_value(M.model)
    return v, M
end

# function optimize(C::Vector, d::Int64, optimizer; kwargs...)
#     X = union([DynamicPolynomials.variables(c[1]) for c in C]...)
#     return optimize(C, X, d, optimizer; kwargs...)
# end
#----------------------------------------------------------------------

