module MOM

import DynamicPolynomials: MonomialVector

using DynamicPolynomials
using MultivariateSeries
using JuMP

function DynamicPolynomials.MonomialVector(V::Vector{PolyVar{true}}, rg::Seq)
    L = DynamicPolynomials.Monomial{true}[]
    for i in rg.val
        append!(L, DynamicPolynomials.monomials(V,i))
    end
    L
end

mutable struct Model
    model::JuMP.Model
end

#----------------------------------------------------------------------
# Define a moment model
function Model(X, d::Int64, optimizer; nu::Int64=1, kwargs...)

    m = JuMP.Model(optimizer;kwargs...)
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
    @variable(m, y[1:s,1:nu])
    m[:moments] = y

    for k in 1:nu 
        H = [ y[m[:index][B[i]*B[j]],k] for i in 1:N, j in 1:N]
        @SDconstraint(m, H >= zeros(N,N))
    end

    return MOM.Model(m)
end

#----------------------------------------------------------------------
function Model(o,  E::Vector, G::Vector, X, d::Int64, optimizer; kwargs...)
    M  = MOM.Model(X,d,optimizer;kwargs...)
    MOM.add_constraint_measure(M)
    MOM.objective(M,o)
    for p in E
        MOM.add_constraint_zero(M,p)
    end
    for p in G
        MOM.add_constraint_nneg(M,p)
    end
    return M
end

#----------------------------------------------------------------------
function Model(pop, X,  d:: Int64, optimizer; kwargs...)
    M = MOM.Model(X, d, optimizer;kwargs...)
    MOM.add_constraint_measure(M)
    for p in constraints(pop)
        if p[2]=="=0"
            MOM.add_constraint_zero(M, p)
        elseif p[2]==">=0"
            MOM.add_constraint_nneg(M, p)
        elseif p[2]=="<=0"
            MOM.add_constraint_nneg(M, -p)
        end
    end
    if objective(pop)[2] == "inf"
        MOM.objective(M,objective(pop)[1])
    else
        MOM.objective(M,-objective(pop)[1])
    end
    return  M
end

#----------------------------------------------------------------------
function add_constraint_zero(M::MOM.Model, eqs...)
    for e in eqs
        p = e*one(Polynomial{true,Float64})
        X = M[:variables]
        L = monomials(X,seq(0:2*M[:degree]-maxdegree(p)))
        for mn in L
            q = p*mn*one(Float64)
            for k in 1:M[:nu]
                @constraint(M.model, sum(t.α*M[:moments][M[:index][t.x],k] for t in q) ==0)
            end
        end
    end
end

#----------------------------------------------------------------------
function add_constraint_nneg(M::MOM.Model, eqs...)
    for e in eqs
        p = e*one(Polynomial{true,Float64})
        X = M[:variables]
        L = monomials(X, seq(0:M[:degree] - maxdegree(p)))
        N = length(L)
        if N == 1
             for k in 1:M[:nu]
                 @constraint(M.model,sum(t.α*M[:moments][M[:index][t.x],k]
                                         for t in p)>=0)
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
end

#----------------------------------------------------------------------
function add_constraint_measure(M::MOM.Model)
    for k in 1:M[:nu]
        @constraint(M.model,M[:moments][1,k] == 1)
    end
end
#----------------------------------------------------------------------
function add_constraint_moments(M::MOM.Model, sigma::Series)
    delta = maxdegree(sigma)
    L = M[:monomials]
    for (m,c) in sigma
        @constraint(M.model,
                    sum(M[:moments][M[:index][m],k]*(-1)^(k-1) for k in 1:M[:nu]) - c==0)
    end
end

#----------------------------------------------------------------------
"""
 Add as objective function the linear functional associated to the polynomial pol to minimize.
"""
function objective(M::MOM.Model, pol)
    f = pol*one(Polynomial{true,Float64})
    @objective(M.model, Min, sum(t.α*M[:moments][M[:index][t.x],1] for t in f))
end

#----------------------------------------------------------------------
function objective_ncl(M)
    B = M[:basis]
    @objective(M.model, Min,  sum(sum(M[:moments][M[:index][B[i]^2],k] for i in 1:length(B)) for k in 1:M[:nu]))
end
#----------------------------------------------------------------------
function objective_tv(M)
    @objective(M.model, Min, sum(M[:moments][1,k] for k in 1:M[:nu]))
end


end  #module MOM
#----------------------------------------------------------------------

#----------------------------------------------------------------------
function getseries(M::MOM.Model)
    [series([M[:monomials][i]=>JuMP.value(M[:moments][i,k])
             for i in 1:length(M[:monomials])])
     for k in 1:M[:nu]]
end

#----------------------------------------------------------------------
function getminimizers(M::MOM.Model)
    s = getseries(M)[1]
    w, Xi = decompose(s);
    Xi
end

#----------------------------------------------------------------------
function getmeasure(M::MOM.Model)
    s = getseries(M)    

    w, Pts = decompose(s[1]);
    for k in 2:M[:nu]
        c, Xi = decompose(s[k])
        w = vcat(w, c*(-1)^(k-1))
        Pts= hcat(Pts,Xi)
    end
    
    return w, Pts
end
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
function JuMP.optimize!(M::MOM.Model)
    JuMP.optimize!(M.model)
end

function JuMP.objective_value(M::MOM.Model)
    JuMP.objective_value(M.model)
end

function JuMP.value(M::MOM.Model)
    JuMP.value.(M[:moments])
end
#----------------------------------------------------------------------
