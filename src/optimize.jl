export optimize, minimize, maximize, mpo_optimizer

#----------------------------------------------------------------------
function optimize(M::JuMP.Model)
    JuMP.optimize!(M)
    v = JuMP.objective_value(M)
    return v, M
end

"""
```julia
v, M = optimize(sense, f, [e1, e2, ...], [p1, p2, ...], X, d)
```
Compute the optimum of `f` under the constraints ``e_i`` =0 and ``p_i \\geq 0`` using a relaxation of order `d` on the moments in the variable `X`. 

  - ``f, e_i, p_i `` should be polynomials in the variables X.
  - 'sense` is a Symbol in [:Inf, :inf, :Min,:min]  or :Sup, :sup, :Max, :max
  - `X` is a tuple of variables
  - `d` is the order of the relaxation
  - `optimizer`is the optimizer used to solve the moment relaxation. By default it is `MPO[:optimizer]`.

If the problem is feasible and has minimizers, it outputs
  - v: the optimum value
  - M: the moment model of type JuMP.Model

Example
-------
```julia
using MomentPolynomialOpt

X  = @polyvar x1 x2
e1 = x1^2-2
e2 = (x2^2-3)*(x1*x2-2)
p1 = x1
p2 = 2-x2
v, M = optimize(:inf, -x1, [e1, e2], [p1, p2], X, 3)
```
To recover the optimizers, see [`get_minimizers`](@ref), [`get_measure`](@ref), [`get_series`](@ref).

"""
function optimize(sense::Symbol, f, Eq::Vector, Pos::Vector,  X, d::Int64, optimizer = MPO[:optimizer]; kwargs...)

    M = MOM.Model(sense, f, Eq, Pos, X, d, optimizer; kwargs...)
    JuMP.optimize!(M)
    v = JuMP.objective_value(M)
    return v, M
end

#------------------------------------------------------------------------
"""
```julia
v, M = minimize(f, [e1, e2, ...], [p1, p2, ...], X, d, optimizer)
```
Compute the infimum of `f` under the constraints ``e_i=0`` and ``p_i \\geq 0`` using a relaxation of order `d` on the moments in the variable `X` and the optimizer `optimizer`.

See [`optimize`](@ref).


"""
function minimize(f, Eq::Vector, Pos::Vector,  X, d::Int64, optimizer=MPO[:optimizer]; kwargs...)
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
function maximize(f, Eq::Vector, Pos::Vector,  X, d::Int64, optimizer=MPO[:optimizer]; kwargs...)
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
using MomentPolynomialOpt, DynamicPolynomials

X  = @polyvar x1 x2
e1 = x1^2-2
e2 = (x2^2-3)*(x1*x2-2)
p1 = x1
p2 = 2-x2
v, M = optimize([(-x1, "inf"), (e1, "=0"), (e2, "=0"), (p1, ">=0"), (p2>=0)], X, 3)
```

To recover the optimal values, see [`get_minimizers`](@ref), [`get_measure`](@ref), [`get_series`](@ref).

"""
function optimize(C::Vector, X, d::Int64, optimizer=MPO[:optimizer]; kwargs...)

    M = MOM.Model(C,X,d,optimizer)
    JuMP.optimize!(M)
    return JuMP.objective_value(M), M
end


#----------------------------------------------------------------------
#import JuMP: set_optimizer

"""
```julia
mpo_optimizer(opt)
```
Define the default optimizer `opt` for the optimization problems created by MomentPolynomialOpt
"""
function mpo_optimizer(opt)
    MPO[:optimizer] = opt 
end


"""
```julia
mpo_optimizer(opt, args...)
```
Define the default optimizer `opt` with its attribute `args...` for the optimization problems created by MomentPolynomialOpt
"""
function mpo_optimizer(opt, args...)
    MPO[:optimizer] = JuMP.optimizer_with_attributes(opt, args...) 
end


#----------------------------------------------------------------------
