module MOM

#import DynamicPolynomials: MonomialVector
using DynamicPolynomials
using MultivariateSeries
using JuMP

mutable struct Model
    model::JuMP.Model
end

    
#----------------------------------------------------------------------
# Define a moment model
function Model(X, d::Int64; nu::Int64=1, w::Vector = [(-1)^(k-1) for k in 1:nu],   kwargs...)

    m = JuMP.Model(kwargs...)

    m[:nu] = nu
    m[:w] = w
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
    @variable(m, y[1:s,1:nu])
    m[:moments] = y

    for k in 1:nu 
        H = [ y[m[:index][B[i]*B[j]],k] for i in 1:N, j in 1:N]
        @SDconstraint(m, H >= zeros(N,N))
    end

    return MOM.Model(m)
end

function Model(X, d::Int64, optimizer; kwargs...)
    M = MOM.Model(X,d; kwargs...)
    set_optimizer(M, optimizer)
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
                        sum(t.α*M[:moments][M[:index][t.x],k] for t in q) == 0)
        end
    end
end
          
#----------------------------------------------------------------------
function add_constraint_zero(M::MOM.Model, idx::AbstractVector, Veq)
    P = [ e*one(Polynomial{true,Float64}) for eq in Veq ]
    X = M[:variables]
    d0 = max([maxdegree(e) for e in P]...)
    L = monomials(X,seq(0:2*M[:degree]-d0))
    for mn in L
        @constraint(M.model,
                    sum(sum(t.α*M[:moments][M[:index][t.x],idx[k]] for t in P[k]*mn) for k in 1:length(idx)) == 0)
    end
end

#----------------------------------------------------------------------
function add_constraint_nneg(M::MOM.Model, pos)
    p = pos*one(Polynomial{true,Float64})
    X = M[:variables]
    d0 = maxdegree(p)
    X = M[:variables]
    L = monomials(X, seq(0:M[:degree] - d0))
    N = length(L)
    if N == 1
        for k in 1:M[:nu]
            @constraint(M.model,
                        sum(t.α*M[:moments][M[:index][t.x],k] for t in p)>=0)
        end
    else
        for k in 1:M[:nu]
            P = [ sum(t.α*M[:moments][M[:index][t.x*L[i]*L[j]],k]
                      for t in p)
                  for i in 1:N, j in 1:N ]
            @SDconstraint(M.model, P >= zeros(N,N))
        end
    end
end

#----------------------------------------------------------------------
function add_constraint_nneg(M::MOM.Model, idx::AbstractVector, Veq)
    p = [ eq*one(Polynomial{true,Float64}) for eq in Veq ]
    X = M[:variables]
    d0 = max([maxdegree(e) for e in p]...)
    L = monomials(X, seq(0:M[:degree] - d0))
    N = length(L)
    if N == 1
        @constraint(M.model,
                    sum(sum(t.α*M[:moments][M[:index][t.x],idx[k]] for t in p[k]) for k in 1:length(idx)) >=0)
    else

        P = [ sum(sum(t.α*M[:moments][M[:index][t.x*L[i]*L[j]],idx[k]]
                      for t in p[k]) for k in 1:length(idx)) 
              for i in 1:N, j in 1:N ]
        @SDconstraint(M.model, P >= zeros(N,N))
    end
end


end  #module MOM

#----------------------------------------------------------------------
function Base.setindex!(p::MOM.Model, v, k::Symbol)  p.model[k] = v end

function Base.getindex(p::MOM.Model, s::Symbol)
    getindex(p.model, s)
end

function Base.show(io::IO, m::MOM.Model)
    println(io, "\nA MOMent program with:")
    Base.show(io, m.model)
end
#----------------------------------------------------------------------
