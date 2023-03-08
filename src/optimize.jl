export optimize, minimize, maximize, mmt_optimizer

#----------------------------------------------------------------------
"""
```julia
v, M = optimize(M)
```
Run the optimizer on the moment program `M` and output the objective_value `v` and the moment program `M`. If the optimization program has no solution, it returns `nothing` and `M`.
"""
function optimize(M::JuMP.Model)
    if haskey(M.obj_dict,:dual)
        JuMP.optimize!(M[:dual])
#        if JuMP.has_values(M[:dual])
            v = JuMP.objective_value(M[:dual])
#        else
#            v = nothing
#        end
    else
        JuMP.optimize!(M)
        v = JuMP.objective_value(M)
    end
    #println("Solver status: ", JuMP.termination_status(M))
    return v, M
end

"""
```julia
v, M = optimize(M, optimizer)
```
Run the optimizer on the moment program `M` using the optimizer `optimizer` and output the objective_value `v` and the moment program `M`. If the optimization program has no solution, it returns `nothing` and `M`. 
"""
function optimize(M::JuMP.Model, optimizer)
    if haskey(M.obj_dict,:dual)
        set_optimizer(M[:dual], optimizer)
        JuMP.optimize!(M[:dual])
#        if JuMP.has_values(M[:dual])
            v = JuMP.objective_value(M[:dual])
#        else
            v = nothing
#        end
    else
        set_optimizer(M, optimizer)
        JuMP.optimize!(M)
        v = JuMP.objective_value(M)
    end
    #println("Solver status: ", JuMP.termination_status(M))
    return v, M
end


#----------------------------------------------------------------------
"""
```julia
v, M = optimize(sense, f, [e1, e2, ...], [p1, p2, ...], X, d)
```
Compute the optimum of `f` under the constraints ``e_i`` =0 and ``p_i \\geq 0`` using a relaxation of order `d` on the moments in the variable `X`. 

  - ``f, e_i, p_i `` should be polynomials in the variables X.
  - 'sense` is a Symbol in [:Inf,:inf,:Min,:min]  or :Sup, :sup, :Max, :max
  - `X` is a tuple of variables
  - `d` is the order of the relaxation
  - `optimizer`is the optimizer used to solve the moment relaxation. By default it is `MMT[:optimizer]`.

If the problem is feasible and has minimizers, it outputs
  - v: the optimum value
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
v, M = optimize(:inf, -x1, [e1, e2], [p1, p2], X, 3)
```
To recover the optimizers, see [`get_minimizers`](@ref), [`get_measure`](@ref), [`get_series`](@ref).

"""
function optimize(sense::Symbol, f, Eq::Vector, Pos::Vector,  X, d::Int64, optimizer = MMT[:optimizer]; kwargs...)

    M = MOM.Model(sense, f, Eq, Pos, X, d, optimizer; kwargs...)
    
    return optimize(M)
end

#------------------------------------------------------------------------
"""
```julia
v, M = minimize(f, [e1, e2, ...], [p1, p2, ...], X, d, optimizer)
```
Compute the infimum of `f` under the constraints ``e_i=0`` and ``p_i \\geq 0`` using a relaxation of order `d` on the moments in the variable `X` and the optimizer `optimizer`.

See [`optimize`](@ref).


"""
function minimize(f, Eq::Vector, Pos::Vector,  X, d::Int64, optimizer=MMT[:optimizer]; kwargs...)
    return optimize(:inf,f,Eq,Pos,X,d,optimizer;kwargs...)
end

#----------------------------------------------------------------------
"""
```julia
v, M = maximize(f, [e1, e2, ...], [p1, p2, ...], X, d, optimizer)
```
Similar to the function `minimize` but compute the supremun of `f`.

See [`optimize`](@ref).
"""
function maximize(f, Eq::Vector, Pos::Vector,  X, d::Int64, optimizer=MMT[:optimizer]; kwargs...)
    return optimize(:sup,f,Eq,Pos,X,d,optimizer;kwargs...)
end

#----------------------------------------------------------------------
"""
```julia
v, M = optimize([(f, set) ...], X, d)
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
function optimize(C::Vector, X, d::Int64, optimizer=MMT[:optimizer]; kwargs...)
    M  = MOM.Model(X,d;kwargs...)
    MOM.constraint_unitmass(M)
    wobj = false
    for c in C
        if c[2] == "inf" || c[2] == "min"
            MOM.set_objective(M, "inf", c[1])
            wobj = true
        elseif c[2] == "sup" || c[2] == "max"
            MOM.set_objective(M, "sup", c[1])
            wobj = true
        elseif c[2] == "=0"
            MOM.constraint_zero(M, c[1])
        elseif c[2] == ">=0"
            MOM.constraint_nneg(M, c[1])
        elseif c[2] == "<=0"
            MOM.constraint_nneg(M,-c[1])
        elseif isa(c[2],AbstractVector)
            MOM.constraint_nneg(M, c[1]-c[2][1])
            MOM.constraint_nneg(M,-c[1]+c[2][2])
        end
    end
    if !wobj
        MOM.set_objective(M, "sup", one(C[1][1]))
    end

    M[:dual] = Dualization.dualize(M,optimizer_with_attributes(optimizer))
    
    return optimize(M)
end


#----------------------------------------------------------------------
#import JuMP: set_optimizer

"""
```julia
mmt_optimizer(opt)
```
Define the default optimizer `opt` for the optimization problems created by MomentTools
"""
function mmt_optimizer(opt)
    MMT[:optimizer] = opt 
end

#----------------------------------------------------------------------
#= """
```julia
v, M = set_optimizer(optimizer)
```
Set the optimizer of the moment program `M` to the dual optimizer of `optimizer`.
"""
function sset_optimizer(M, optimizer)
    if haskey(M,:type) && M[:type] == :moment
        @info "Using dual optimizer for the model"
        JuMP.set_optimizer(M, Dualization.dual_optimizer(optimizer))
    else
        set_optimizer(M, optimizer)
    end
end
=#
