module MOS

import DynamicPolynomials: MonomialVector

using DynamicPolynomials
using MultivariateSeries
using JuMP

mutable struct Model
    model::JuMP.Model
end

# Define a moment model
function Model(X, d::Int64, optimizer; kwargs...)
    m = JuMP.Model(optimizer;kwargs...)
    m[:variables] = X

    m[:degree] = d
    B = monomials(X,seq(0:d))
    N = length(B)

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
    @variable(m, y[1:s])
    m[:moments] = y

    @constraint(m, m[:moments][1]==1)
    
    m[:H] = [ y[m[:index][B[i]*B[j]]] for i in 1:N, j in 1:N]
    @SDconstraint(m, m[:H] >= zeros(N,N))

    return MOS.Model(m)
end

#----------------------------------------------------------------------
function Model(o,  E::Vector, G::Vector, X, d::Int64, optimizer; kwargs...)
    m  = MOS.Model(X,d,optimizer;kwargs...)
    MOS.objective(m,o)
    for p in E
        MOS.add_constraint_zero(m,p)
    end
    for p in G
        MOS.add_constraint_nneg(m,p)
    end
    return m
end


#----------------------------------------------------------------------

function Model(pop, X,  d:: Int64, optimizer; kwargs...)
    m = MOS.Model(X, d, optimizer)
    for p in constraints(pop)
        if p[2]=="=0"
            MOS.add_constraint_zero(m, p)
        elseif p[2]==">=0"
            MOS.add_constraint_nneg(m, p)
        elseif p[2]=="<=0"
            MOS.add_constraint_nneg(m, -p)
        end
    end
    if objective(pop)[2] == "inf"
        MOS.objective(m,objective(pop)[1])
    else
        MOS.objective(m,-objective(pop)[1])
    end
    return  m
end

#----------------------------------------------------------------------
function add_constraint_zero(m::MOS.Model, eqs...)
    for e in eqs
        p = e*one(Polynomial{true,Float64})
        X = m.model[:variables]
        mon = monomials(X,seq(0:2*m.model[:degree]-maxdegree(p)))
        for mn in mon
            q = p*mn*one(Float64)
            @constraint(m.model, sum(t.α*m.model[:moments][m.model[:index][t.x]] for t in q) ==0)
        end
    end
end

#----------------------------------------------------------------------
function add_constraint_nneg(m::MOS.Model, eqs...)
    for e in eqs
        p = e*one(Polynomial{true,Float64})
        X = m.model[:variables]
        L = monomials(X, seq(0:m[:degree] - maxdegree(p)))
        N = length(L)
        P = [ sum(t.α*m.model[:moments][m.model[:index][t.x*L[i]*L[j]]] for t in p)
              for i in 1:N, j in 1:N ]
        @SDconstraint(m.model, P >= zeros(N,N))
    end
end

#----------------------------------------------------------------------
function add_constraint_moments(m::MOS.Model, sigma::Series)
    delta = maxdegree(sigma)
    L = m[:monomials]
    N = length(L)
    for i in 1:N
        for j in 1:N
            mn = L[i]*L[j]
	    cf = sigma[mn]
	    if  maxdegree(mn) <= delta
                @constraint(m, sum(H[i,j]*s for (H,s) in zip(m[:H], m[:sign]))-cf == 0)
            end
        end
    end
end


#----------------------------------------------------------------------
"""
 Add as objective function the linear functional associated to the polynomial pol
"""
function objective(m::MOS.Model, pol)
    f = pol*one(Polynomial{true,Float64})

    @objective(m.model, Min, sum(t.α*m.model[:moments][m.model[:index][t.x]] for t in f))
end

end

#----------------------------------------------------------------------


function getseries(m::MOS.Model)
    v = JuMP.value.(m.model[:moments])
    s = series([m.model[:monomials][i]=> v[i] for i in 1:length(v)])
end

#----------------------------------------------------------------------
function getminimizers(m::MOS.Model)
    s = getseries(m)
    w, Xi = decompose(s);
    Xi
end

function Base.setindex!(p::MOS.Model, v, k::Symbol)  p.model[k] = v end

function Base.getindex(p::MOS.Model, s::Symbol)
    getindex(p.model, s)
end

function Base.show(io::IO, m::MOS.Model)
    println(io, "\nA MOment sequenceS program with:")
    Base.show(io, m.model)
end
