module MOM

#import DynamicPolynomials: MonomialVector
using DynamicPolynomials
using MultivariateSeries
using JuMP
using LinearAlgebra

mutable struct Model
    model::JuMP.Model
end

    
#----------------------------------------------------------------------
# Define a moment model
function Model(X, d::Int64; nu::Int64=1,   kwargs...)

    M = JuMP.Model(kwargs...)

    M[:nu] = nu
    M[:variables] = X
    M[:degree] = d
    
    B = monomials(X,seq(0:d))
    N = length(B)
    M[:basis] = B

    M[:index] = Dict{DynamicPolynomials.Monomial{true},Tuple{Int64,Int64}}()

    M[:monomials] = DynamicPolynomials.Monomial{true}[]

    M[:H] = []
    for k in 1:nu
        push!(M[:H], @variable(M, [1:N,1:N], SDP))
    end

    # Hankel structure
    for i in 1:N
        for j in 1:N
            mn = M[:basis][i]*M[:basis][j]
	        if !haskey(M[:index], mn)
	            M[:index][mn] = (i,j)
                    push!(M[:monomials], mn)
	        else
                    id = get(M[:index], mn, (0,0))
	            if (i != id[2]) || (j != id[1])
                        for k in 1:nu
	                    @constraint(m, M[:H][k][i,j]-M[:H][k][id...] == 0)
                        end
                    end
	        end
	end
    end

    return MOM.Model(M)
end

function Model(X, d::Int64, optimizer; kwargs...)
    M = MOM.Model(X,d; kwargs...)
    set_optimizer(M.model, optimizer)
    return M
end


#----------------------------------------------------------------------
function add_constraint_zero(M::MOM.Model, eq)
    p = eq*one(Polynomial{true,Float64})
    X = M[:variables]
    L = monomials(X,seq(0:2*M[:degree]-maxdegree(p)))

    for mn in L
        for k in 1:M[:nu]
            @constraint(M.model,
                        sum(t.α*M[:H][k][M[:index][t.x]...] for t in p*mn)==0)
        end
    end
end
          
#----------------------------------------------------------------------
function add_constraint_zero(M::MOM.Model, idx::Vector{Int64}, Veq)
    P = [ e*one(Polynomial{true,Float64}) for eq in Veq ]
    X = M[:variables]
    d0 = max([maxdegree(e) for e in P]...)
    L = monomials(X,seq(0:2*M[:degree]-d0))
    for mn in L
        @constraint(M.model,
                    sum(sum(t.α*M[:H][idx[k]][M[:index][t.x]...] for t in P[k]*mn) for k in 1:length(idx)) == 0)
    end
end

#----------------------------------------------------------------------
function add_constraint_nneg(M::MOM.Model, pos)
    p = pos*one(Polynomial{true,Float64})
    X = M[:variables]
    d0 = Int64(ceil(maxdegree(p)/2))
    X = M[:variables]
    L = monomials(X, seq(0:M[:degree] - d0))
    N = length(L)
    if N == 1
        for k in 1:M[:nu]
            @constraint(M.model,
                        sum(t.α*M[:moments][M[:index][t.x],k] for t in p)>=0)
        end
    else
        P  = @variable(m, [1:N,1:N], SDP)
        for i in 1:N
            for j  in 1:i
                mn = L[i]*L[j]
                for k in 1:M[:nu]
                    @constraint(M.model,
                                P[i,j]-sum(t.α*M[:H][k][M[:index][t.x*mn]...] for t in p) == 0)
                end
            end
        end
    end
end

#----------------------------------------------------------------------
function add_constraint_nneg(M::MOM.Model, idx::Vector{Int64}, Veq)
    p = [ eq*one(Polynomial{true,Float64}) for eq in Veq ]
    X = M[:variables]
    d0 = Int64(max([maxdegree(e) for e in p]...)/2)
    L = monomials(X, seq(0:M[:degree] - d0))
    N = length(L)
    if N == 1
        @constraint(M.model,
                    sum(sum(t.α*M[:H][idx[k]][M[:index][t.x]...] for t in p[k]) for k in 1:length(idx)) >=0)
    else
        P  = @variable(m, [1:N,1:N], SDP)
        for i in 1:N
            for j  in 1:i
                mn = L[i]*L[j]
                @constraint(M.model,
                            P[i,j]-sum(sum(t.α*M[:H][idx[k]][M[:index][t.x*L[i]*L[j]]...] for t in p[k]) for k in 1:length(idx)) == 0 ) 
            end
        end
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
