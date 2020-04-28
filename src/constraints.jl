using DynamicPolynomials
export constraint_zero, constraint_nneg, constraint_unitmass, constraint_moments,
    objective, objective_tv, objective_ncl


#----------------------------------------------------------------------
function constraint_zero(M::MOM.Model, idx::Vector{Int64}, eqs::Vector...)
    for e in eqs
        MOM.add_constraint_zero(M,idx,e)
    end
end


function constraint_zero(M::MOM.Model, eqs...)
    for e in eqs
        MOM.add_constraint_zero(M,e)
    end
end

#----------------------------------------------------------------------
function constraint_nneg(M::MOM.Model, idx::Vector{Int64}, eqs::Vector...)
    for e in eqs
        MOM.add_constraint_nneg(M,idx,e)
    end
end

function constraint_nneg(M::MOM.Model, eqs...)
    for e in eqs
        MOM.add_constraint_nneg(M,e)
    end
end

#----------------------------------------------------------------------
function constraint_moments(M::MOM.Model, moments::Vector, P::Vector=[(-1)^(k-1) for k in 1:M[:nu]])
    for c in moments
        p = [ q*c[1]*one(Polynomial{true,Float64}) for q in P ]
        @constraint(M.model,
                    sum(sum(t.α*M[:moments][M[:index][t.x],k] for t in p[k]) for k in 1:M[:nu])-c[2]==0)
    end
end

#----------------------------------------------------------------------
function constraint_moments(M::MOM.Model, moments::Vector, idx::Vector{Int64}, P::Vector)
    for c in moments
        p = [q*c[1]*one(Polynomial{true,Float64}) for q in P]
        @constraint(M.model,
                    sum(sum(t.α*M[:moments][M[:index][t.x],idx[k]] for t in p[k]) for k in 1:length(idx))-c[2]==0)
    end
end


#----------------------------------------------------------------------
function constraint_unitmass(M::MOM.Model, P::Vector=[(-1)^(k-1) for k in 1:M[:nu]])
    p = [ q*one(Polynomial{true,Float64}) for q in P ]
    @constraint(M.model,
                sum(sum(t.α*M[:moments][M[:index][t.x],k] for t in p[k]) for k in 1:M[:nu])-1==0)
end

#----------------------------------------------------------------------
function constraint_unitmass(M::MOM.Model, idx::Vector{Int64}, P::Vector)
    p = [q*one(Polynomial{true,Float64}) for q in P]
    @constraint(M.model,
                sum(sum(t.α*M[:moments][M[:index][t.x],idx[k]] for t in p[k]) for k in 1:length(idx))-1==0)
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
"""
 Add as objective function the linear functional associated to the polynomial pol to minimize.
"""
function objective(M::MOM.Model, idx::Vector{Int64}, pol::Vector, sense="inf")
    if pol != nothing
        f = [p*one(Polynomial{true,Float64}) for p in pol]
        if sense == "inf"
            @objective(M.model, Min, sum(sum(t.α*M[:moments][M[:index][t.x],k]*M[:w][idx[k]] for t in f[k]) for k in 1:length(idx)))
        else
            @objective(M.model, Max, sum(sum(t.α*M[:moments][M[:index][t.x],k]*M[:w][idx[k]] for t in f[k]) for k in 1:length(idx)))
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
