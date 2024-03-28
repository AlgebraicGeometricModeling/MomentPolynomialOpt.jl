export MOM
import JuMP: optimize!, objective_value

module MOM


using DynamicPolynomials
#using MultivariateSeries
using JuMP, Dualization

# import MathOptInterface
# const MOI = MathOptInterface

using LinearAlgebra


import MomentPolynomialOpt:MMT

convert_Float64 = function(pol)
    if typeof(pol) != Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Int64}
        return dot(Float64.(coefficients(pol)),monomials(pol))
    else
        return pol
    end
end

export Model
#----------------------------------------------------------------------
"""
Construct the Moment Program in the variables X of order d.
The moments of all monomials in X of degree 2*d are variables of
the optimization program.

```
M = MOM.Model(X,d, optimizer=MMT[:optimizer]; nu=k)
```
  - `X` is the vector of variables
  - `d` is the order of the moment relaxation
  - 'optimizer`is the the optimizer used to solve the SDP program
  - `nu=k` is the number of Positive Moment Sequences
"""
function Model(X, d::Int64, optimizer=MMT[:optimizer]; nu::Int64=1, kwargs...)

    M = JuMP.Model(kwargs...)

    M[:type] = :moment
    M[:nu] = nu
    M[:variables] = X
    M[:degree] = d


    B = monomials(X,0:d)
    N = length(B)
    M[:basis] = B

    M[:monomials] = typeof(B[1])[]
    M[:index] = Dict{typeof(B[1]),Int64}()

    # Hankel structure
    c = 0
    for i in 1:N
        for j  in 1:N
            mn = B[i]*B[j]
	    if !haskey(M[:index], mn)
                c += 1
                M[:index][mn] = c
                push!(M[:monomials], mn)
	    end
        end
    end

    s = length(M[:monomials])
    @variable(M, y[1:nu, 1:s])
    M[:moments] = y

    for k in 1:nu
        H = [ y[k,M[:index][B[i]*B[j]]]+0 for i in 1:N, j in 1:N]
        @constraint(M, Symmetric(H) in PSDCone())
    end
    
#  JuMP.set_optimizer(M, JuMP.optimizer_with_attributes(optimizer))
    #JuMP.set_optimizer(M, JuMP.optimizer_with_attributes(DualOptimizer,optimizer()))
    #@info "Using dual optimizer"
    
    return M #MOM.Model(m)
end

#----------------------------------------------------------------------

function dualize!(M::JuMP.Model, optimizer=MMT[:optimizer])
    if optimizer == nothing
        M[:dual] = Dualization.dualize(M)
    else
        M[:dual] = Dualization.dualize(M,optimizer); #optimizer_with_attributes(optimizer))
    end
end

#----------------------------------------------------------------------

include("constraints.jl")
include("objective.jl")

#----------------------------------------------------------------------
"""
Construct the Moment Program in the variables X of order d.
The moments of all monomials in X of degree 2*d are variables of
the optimization program.

```
M = MomentModel(X,d; nu=k)
```
  - `X` is the vector of variables
  - `d` is the order of the moment relaxation.
  - `nu=k` is the number of Positive Moment Sequences
"""
function MomentModel(X, d::Int64; nu::Int64=1,  kwargs...)
    return MOM.Model(X,d;nu=nu,kwargs...)
end

function MomentModel(X, d::Int64, optimizer; kwargs...)
    return MOM.Model(X,d,optimizer; kwargs...)
end


"""
```julia
M = MOM.Model( `sense`, f, [e1, e2, ...], [g1, g2, ...], X, d)
```
Construct the Moment Program in the variables X of order d.
   - `sense` == "inf" or "sup"
   - `f` polynomial objective function
   - `[e1, e2, ...]` array of polynomial equality constraints (can be empty)
   - `[g1, g2, ...]` array of non-negativity constraints (can be empty)
   - `X` is the vector of variables
   - `d` is the order of the moment relaxation.
"""
function Model(sense::Symbol, f, Eq::Vector, Pos::Vector,  X, d::Int64, optimizer=MMT[:optimizer])
    M = MOM.Model(X, d, optimizer)
    constraint_unitmass(M)
    for e in Eq  constraint_zero(M, e) end
    for p in Pos constraint_nneg(M, p) end
    if f != nothing
        if in(sense,[:Inf,:inf,:Min,:min])
            set_objective(M, "inf", f)
        else
            set_objective(M, "sup", f)
        end
    else
        set_objective(M, "sup", one(Polynomial{true,Float64}))
    end

    MOM.dualize!(M,optimizer)
     
    return M
end



"""
```julia
M = MOM.Model(C, X, d)
```
Construct the Moment Program where
   - C is a vector of pairs (f, sense ) of objective or constraints where f is a polynomial and sense is "inf", "min", "sup", "max", ">=0", "<=0", "=0", or an interval 
   - `X` is the vector of variables
   - `d` is the order of the moment relaxation.
"""
function  Model(C::Vector, X, d::Int64, optimizer = MMT["optimizer"]; kwargs...)
    M = MOM.Model(X, d, optimizer; kwargs...)
    constraint_unitmass(M)
    for c in C
        if c[2] == "inf" || c[2] == "min"
            set_objective(M, "inf", c[1])
            wobj = true
        elseif c[2] == "sup" || c[2] == "max"
            set_objective(M,"sup", c[1])
            wobj = true
        elseif c[2] == "=0"
            constraint_zero(M, c[1])
        elseif c[2] == ">=0"
            constraint_nneg(M, c[1])
        elseif c[2] == "<=0"
            constraint_nneg(M, -c[1])
        elseif isa(c[2], AbstractVector)
            constraint_nneg(M, c[1] - c[2][1])
            constraint_nneg(M, -c[1] + c[2][2])
        end
    end

    MOM.dualize!(M, optimizer)
    return M
end

end  #module MOM

