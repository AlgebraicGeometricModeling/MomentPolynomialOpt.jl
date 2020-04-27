using DynamicPolynomials
export constraint_zero, constraint_nneg, constraint_unitmass, constraint_moments,
    objective, objective_tv, objective_ncl


#----------------------------------------------------------------------
function constraint_zero(M::MOM.Model, idx::Vector{Int64}, eqs::Polynomial...)
    for e in eqs
        MOM.add_constraint_zero(M,idx,e)
    end
end


function constraint_zero(M::MOM.Model, eqs::Polynomial...)
    constraint_zero(M,collect(1:M[:nu]),eqs...)
end

#----------------------------------------------------------------------
function constraint_nneg(M::MOM.Model, idx::Vector{Int64}, eqs::Polynomial...)
    for e in eqs
        MOM.add_constraint_nneg(M,idx,e)
    end
end

function constraint_nneg(M::MOM.Model, eqs::Polynomial...)
    constraint_nneg(M,collect(1:M[:nu]),eqs...)
end
#----------------------------------------------------------------------
function constraint_unitmass(M::MOM.Model, idx=1:M[:nu])
    for k in idx
        @constraint(M.model,M[:moments][1,k] - 1 == 0)
    end
end

function constraint_unitmass(M::MOM.Model, p::Polynomial)
    for k in 1:M[:nu]
         @constraint(M.model,
                    sum(t.α*M[:moments][M[:index][t.x],k] for t in p) -1 ==0)
    end
end

#----------------------------------------------------------------------
function constraint_moments(M::MOM.Model, sigma::Series)
    delta = maxdegree(sigma)
    L = M[:monomials]
    for (m,c) in sigma
        @constraint(M.model,
                    sum(M[:moments][M[:index][m],k]*M[:w][k] for k in 1:M[:nu]) - c == 0)
    end
end

#----------------------------------------------------------------------
function constraint_moments(M::MOM.Model, moments::Vector)
    for c in moments
        p = c[1]*one(Polynomial{true,Float64}) 
        @constraint(M.model,
                    sum(sum(t.α*M[:moments][M[:index][t.x],k] for t in p)*M[:w][k] for k in 1:M[:nu])-c[2]==0)
    end
end

#----------------------------------------------------------------------
function constraint_moments(M::MOM.Model, moments::Vector, k::Int64)
    for c in moments
        p = c[1]*one(Polynomial{true,Float64}) 
        @constraint(M.model,
                    sum(t.α*M[:moments][M[:index][t.x],k] for t in p)*M[:w][k]-c[2]==0)
    end
end

#----------------------------------------------------------------------
"""
 Add as objective function the linear functional associated to the polynomial pol to minimize.
"""
function objective(M::MOM.Model, pol, sense="inf")
    if pol != nothing
        f = pol*one(Polynomial{true,Float64})
        if sense == "inf"
            @objective(M.model, Min, sum(sum(t.α*M[:moments][M[:index][t.x],k]*M[:w][k] for t in f) for k in 1:M[:nu]))
        else
            @objective(M.model, Max, sum(sum(t.α*M[:moments][M[:index][t.x],k]*M[:w][k] for t in f) for k in 1:M[:nu]))
        end
    end
end

#----------------------------------------------------------------------
function objective_ncl(M)
    B = M[:basis]
    @objective(M.model, Min,  sum(sum(M[:moments][M[:index][B[i]^2],k] for i in 1:length(B)) for k in 1:M[:nu]))
end
#----------------------------------------------------------------------
function objective_tv(M)
    @objective(M.model, Min, sum(M[:moments][1,k]*abs(M[:w][k]) for k in 1:M[:nu]))
end
