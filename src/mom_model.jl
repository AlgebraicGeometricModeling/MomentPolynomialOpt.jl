export MOM
import JuMP: optimize!, objective_value

module MOM


using DynamicPolynomials
using MultivariateSeries
using JuMP

# import MathOptInterface
# const MOI = MathOptInterface


using LinearAlgebra
import Dualization

#mutable struct Model
#    model::JuMP.Model
#end


import MomentTools:MMT

#----------------------------------------------------------------------
"""
Construct the Moment Program in the variables X of order d.
The moments of all monomials in X of degree 2*d are variables of
the optimization program.

```
M = MOM.Model(X,d; nu=k)
```
  - `X` is the vector of variables
  - `d` is the order of the moment relaxation.
  - `nu=k` is the number of Positive Moment Sequences
"""
function Model(X, d::Int64; nu::Int64=1, dual::Bool=true,  kwargs...)

    M = JuMP.Model(kwargs...)

    M[:dual] = dual
    M[:nu] = nu
    M[:variables] = X
    M[:degree] = d

    if haskey(MMT,:optimizer,)
        JuMP.set_optimizer(M, Dualization.dual_optimizer(MMT[:optimizer]))
    end

    B = monomials(X,seq(0:d))
    N = length(B)
    M[:basis] = B

    M[:monomials] = DynamicPolynomials.Monomial{true}[]
    M[:index] = Dict{DynamicPolynomials.Monomial{true},Int64}()

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
        #@SDconstraint(m, H >= zeros(N,N))
    end

    return M #MOM.Model(m)
end


function Model(X, d::Int64, optimizer; kwargs...)
    M = MOM.Model(X,d; kwargs...)
    if M[:dual]
        set_optimizer(M, Dualization.dual_optimizer(optimizer))
    else
        set_optimizer(M, optimizer)
    end
    return M
end


#----------------------------------------------------------------------
#function get_series(M::JuMP.Model)
#    [series([M[:monomials][i]=>JuMP.value(M[:moments][k,i])
#             for i in 1:length(M[:monomials])]) for k in 1:M[:nu]]
#end


include("MOM/constraints.jl")
include("MOM/objective.jl")
include("MOM/optimize.jl")



#----------------------------------------------------------------------
#= function Base.setindex!(p::JuMP.Model, v, k::Symbol)  p.model[k] = v end

function Base.getindex(p::JuMP.Model, s::Symbol)
    getindex(p.model, s)
end

function Base.show(io::IO, m::JuMP.Model)
    println(io, "\nA Moment program with:")
    Base.show(io, m.model)
end
=#
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
function Model(sense::Symbol, f, Eq::Vector, Pos::Vector,  X, d::Int64; kwargs...)
    M = MOM.Model(X, d; kwargs...)
    constraint_unitmass(M)
    constraint_zero(M,Eq...)
    constraint_nneg(M,Pos...)
    if f != nothing
        if sense == :Inf
            set_objective(M, "inf", f)
        else
            set_objective(M, "sup", f)
        end
    else
        set_objective(M, "sup", one(Polynomial{true,Float64}))
    end
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
function  Model(C::Vector, X, d::Int64; kwargs...)
    M = MOM.Model(X, d; kwargs...)
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
    return M
end

#function JuMP.optimize!(M::MOM.Model)
#    JuMP.optimize!(M)
#    return JuMP.objective_value(M)
#end

end  #module MOM