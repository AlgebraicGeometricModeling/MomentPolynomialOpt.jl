
function optimize(sense::Symbol, fct, Eq::Vector, Pos::Vector,  X, d::Int64; kwargs...)
    M = MOM.Model(X, d; kwargs...)
    constraint_unitmass(M)
    constraint_zero(M,Eq...)
    constraint_nneg(M,Pos...)
    if fct != nothing
       if sense == :Inf
           set_objective(M, "inf", fct)
       else
           set_objective(M, "sup", fct)
       end
    else
        set_objective(M, "inf", one(Polynomial{true,Float64}))
    end
    JuMP.optimize!(M)
    v = JuMP.objective_value(M)
    return v, M
end

#----------------------------------------------------------------------
"""
```julia
v, M = MOM.set_optimizer(optimizer)
```
Set the optimizer of the moment program `M` to the dual optimizer of `optimizer`.
"""
function set_optimizer(optimizer)

#    if M[:dual]
        MOMENV[:optimizer] = Dualization.dual_optimizer(optimizer)
#    else
#        set_optimizer(M, optimizer)
#    end
end

function set_optimizer(M, optimizer)
    if M[:dual]
        JuMP.set_optimizer(M, Dualization.dual_optimizer(optimizer))
    else
        set_optimizer(M, optimizer)
    end
end
#------------------------------------------------------------------------
"""
```julia
v, M = MOM.minimize(f, [e1, e2, ...], [p1, p2, ...], X, d, optimizer)
```
Compute the infimum of `f` under the constraints ``e_i`` =0 and ``p_i \\geq 0`` using a relaxation of order `d` on the moments in the variable `X` and the optimizer `M`.

``f, e_i, p_i `` should be polynomials in the variables X.

If the problem is feasible and has minimizers, it outputs
  - v: the minimum value
  - M: the moment model of type JuMP.Model

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
        set_objective(M, "inf", fct)
    else
        set_objective(M, "inf", one(Polynomial{true,Float64}))
    end
    JuMP.optimize!(M)
    v = objective_value(M)
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
        set_objective(M, "sup", fct)
    else
        set_objective(M, "sup", one(Polynomial{true,Float64}))
    end
    JuMP.optimize!(M)
    v = objective_value(M)
    return v, M
end

#----------------------------------------------------------------------
"""
```julia
v, M = optimize([(f, set), ...], X, d)
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
function optimize(C::Vector, X, d::Int64; kwargs...)
    M  = MOM.Model(X,d;kwargs...)
    constraint_unitmass(M)
    wobj = false
    for c in C
        if c[2] == "inf" || c[2] == "min"
            set_objective(M, "inf", c[1])
            wobj = true
        elseif c[2] == "sup" || c[2] == "max"
            set_objective(M, "sup", c[1])
            wobj = true
        elseif c[2] == "=0"
            constraint_zero(M, c[1])
        elseif c[2] == ">=0"
            constraint_nneg(M, c[1])
        elseif c[2] == "<=0"
            constraint_nneg(M,-c[1])
        elseif isa(c[2],AbstractVector)
            constraint_nneg(M, c[1]-c[2][1])
            constraint_nneg(M,-c[1]+c[2][2])
        end
    end
    if !wobj
        set_objective(M, "sup", one(C[1][1]))
    end
    JuMP.optimize!(M)
    v = objective_value(M)
    return v, M
end

