export MOM
import JuMP: optimize!, objective_value

module MOM

#import DynamicPolynomials: MonomialVector
using DynamicPolynomials
using MultivariateSeries
using JuMP

# import MathOptInterface
# const MOI = MathOptInterface


using LinearAlgebra
import Dualization

mutable struct Model
    model::JuMP.Model
end


#----------------------------------------------------------------------
# Define a moment model
function Model(X, d::Int64; nu::Int64=1, dual::Bool=true,  kwargs...)

    m = JuMP.Model(kwargs...)

    m[:dual] = dual
    m[:nu] = nu
    m[:variables] = X
    m[:degree] = d

    B = monomials(X,seq(0:d))
    N = length(B)
    m[:basis] = B

    m[:monomials] = DynamicPolynomials.Monomial{true}[]
    m[:index] = Dict{DynamicPolynomials.Monomial{true},Int64}()

    # Hankel structure
    c = 0
    for i in 1:N
        for j  in 1:N
            mn = B[i]*B[j]
	    if !haskey(m[:index], mn)
                c += 1
                m[:index][mn] = c
                push!(m[:monomials], mn)
	    end
        end
    end

    s = length(m[:monomials])
    @variable(m, y[1:nu, 1:s])
    m[:moments] = y

    for k in 1:nu
        H = [ y[k,m[:index][B[i]*B[j]]]+0 for i in 1:N, j in 1:N]
        @constraint(m, Symmetric(H) in PSDCone())
        #@SDconstraint(m, H >= zeros(N,N))
    end

    return MOM.Model(m)
end


function Model(X, d::Int64, optimizer; kwargs...)
    M = MOM.Model(X,d; kwargs...)
    if M[:dual]
        set_optimizer(M.model, Dualization.dual_optimizer(optimizer))
    else
        set_optimizer(M.model, optimizer)
    end
    return M
end

#----------------------------------------------------------------------
function add_constraint_zero(M::MOM.Model, eq)
    p = eq*one(Polynomial{true,Float64})
    X = M[:variables]
    L = monomials(X,seq(0:2*M[:degree]-maxdegree(p)))
    for mn in L
        q = p*mn*one(Float64)
        for k in 1:M[:nu]
            @constraint(M.model,
                        sum(t.α*M[:moments][k,M[:index][t.x]] for t in q) == 0)
        end
    end
end

function add_constraint_zero(M::MOM.Model, idx::Vector{Int64}, Veq)
    P = [ e*one(Polynomial{true,Float64}) for eq in Veq ]
    X = M[:variables]
    d0 = max([maxdegree(e) for e in P]...)
    L = monomials(X,seq(0:2*M[:degree]-d0))
    for mn in L
        @constraint(M.model,
                    sum(sum(t.α*M[:moments][idx[k],M[:index][t.x]] for t in P[k]*mn) for k in 1:length(idx)) == 0)
    end
end

#----------------------------------------------------------------------
function add_constraint_nneg(M::MOM.Model, pos)
    p = pos*one(Polynomial{true,Float64})
    d0 = Int(ceil(maxdegree(p)/2))
    X = M[:variables]
    L = monomials(X, seq(0:M[:degree] - d0))
    N = length(L)
    if N == 1
        for k in 1:M[:nu]
            @constraint(M.model,
                        sum(t.α*M[:moments][k,M[:index][t.x]] for t in p)>=0)
        end
    else
        for k in 1:M[:nu]
            P = [ sum(t.α*M[:moments][k,M[:index][t.x*L[i]*L[j]]]
                      for t in p)+0
                  for i in 1:N, j in 1:N ]
            @constraint(M.model, Symmetric(P) in PSDCone())
        end
    end
end

function add_constraint_nneg(M::MOM.Model, idx::Vector{Int64}, Veq)
    p = [ eq*one(Polynomial{true,Float64}) for eq in Veq ]
    X = M[:variables]
    d0 = Int(ceil(max([maxdegree(e) for e in p]...)/2))
    L = monomials(X, seq(0:M[:degree] - d0))
    N = length(L)
    if N == 1
        @constraint(M.model,
                    sum(sum(t.α*M[:moments][idx[k],M[:index][t.x]] for t in p[k]) for k in 1:length(idx)) >=0)
    else
        P = [ sum(sum(t.α*M[:moments][idx[k],M[:index][t.x*L[i]*L[j]]]
                      for t in p[k]) for k in 1:length(idx))+0
              for i in 1:N, j in 1:N ]
        @constraint(M.model, Symmetric(P) in PSDCone())
    end
end

#----------------------------------------------------------------------
function add_constraint_moment(M::MOM.Model, v , p)
    for k in 1:M[:nu]
        @constraint(M.model,
                    sum(sum(t.α*M[:moments][k,M[:index][t.x]] for t in p) for k in 1:M[:nu])-v == 0)
    end
end

function add_constraint_moment(M::MOM.Model, v, p::Vector)
    @constraint(M.model,
                sum(sum(t.α*M[:moments][k,M[:index][t.x]] for t in p[k]) for k in 1:length(p)) - v ==0)
end

function add_constraint_moment(M::MOM.Model, v, idx::Vector{Int64}, p::Vector)
    @constraint(M.model,
                sum(sum(t.α*M[:moments][idx[k],M[:index][t.x]] for t in p[k]) for k in 1:length(idx)) - v ==0)
end


#----------------------------------------------------------------------
function get_series(M::MOM.Model)
    [series([M[:monomials][i]=>JuMP.value(M[:moments][k,i])
             for i in 1:length(M[:monomials])]) for k in 1:M[:nu]]
end

end  #module MOM
#----------------------------------------------------------------------

#----------------------------------------------------------------------
function Base.setindex!(p::MOM.Model, v, k::Symbol)  p.model[k] = v end

function Base.getindex(p::MOM.Model, s::Symbol)
    getindex(p.model, s)
end

function Base.show(io::IO, m::MOM.Model)
    println(io, "\nA Moment program with:")
    Base.show(io, m.model)
end
#----------------------------------------------------------------------
export MomentModel

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
M = MomentModel((`sense`,f), [e1, e2, ...], [g1, g2, ...], X, d)
```
Construct the Moment Program in the variables X of order d.
   - `sense` == "inf" or "sup"
   - `f` polynomial objective function
   - `[e1, e2, ...]` array of polynomial equality constraints (can be empty)
   - `[g1, g2, ...]` array of non-negativity constraints (can be empty)
   - `X` is the vector of variables
   - `d` is the order of the moment relaxation.
"""
function MomentModel(fct, Eq::Vector, Pos::Vector,  X, d::Int64; kwargs...)
    M = MOM.Model(X, d; kwargs...)
    constraint_unitmass(M)
    constraint_zero(M,Eq...)
    constraint_nneg(M,Pos...)
    if fct != nothing
        set_objective(M, fct[1], fct[2])
    else
        set_objective(M, "sup", one(Polynomial{true,Float64}))
    end
    return M
end

"""
```julia
M = MomentModel(C, X, d)
```
Construct the Moment Program where
   - C is a vector of pairs (f, sense ) of objective or constraints where f is a polynomial and sense is "inf", "min", "sup", "max", ">=0", "<=0", "=0", or an interval 
   - `X` is the vector of variables
   - `d` is the order of the moment relaxation.
"""
function  MomentModel(C::Vector, X, d::Int64; kwargs...)
    M = MOM.Model(X, d; kwargs...)
    constraint_unitmass(M)
    for c in C
        if c[2] == "inf" || c[2] == "min"
            MomentTools.objective(M, c[1], "inf")
            wobj = true
        elseif c[2] == "sup" || c[2] == "max"
            MomentTools.objective(M, c[1], "sup")
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

function JuMP.optimize!(M::MOM.Model)
    JuMP.optimize!(M.model)
    return JuMP.objective_value(M.model)
end

