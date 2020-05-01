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
function constraint_moments(M::MOM.Model, moments::Vector)
    for c in moments
        MOM.add_constraint_moment(M, c[2], c[1]*one(Polynomial{true,Float64}))
    end
end

#----------------------------------------------------------------------
function constraint_moments(M::MOM.Model, moments::Vector, idx::Vector{Int64}, P::Vector)
    for c in moments
        MOM.add_constraint_moment(M, c[2], idx, P*c[1]*one(Polynomial{true,Float64}))
    end
end


#----------------------------------------------------------------------
function constraint_unitmass(M::MOM.Model)
    MOM.add_constraint_moment(M, 1, one(Polynomial{true,Float64}))
end

#----------------------------------------------------------------------
function constraint_unitmass(M::MOM.Model, idx::Vector{Int64}, P::Vector)
    MOM.add_constraint_moment(M, 1, idx, P)
end

#----------------------------------------------------------------------
"""
 Add as objective function the linear functional associated to the polynomial pol to minimize.
"""
function objective(M::MOM.Model, pol, sense="inf")
    if pol != nothing
        f = pol*one(Polynomial{true,Float64})
        MOM.set_objective(M,f,sense)
    end
end

#----------------------------------------------------------------------
"""
 Add as objective function the linear functional associated to the polynomial pol to minimize.
"""
function objective(M::MOM.Model, idx::Vector{Int64}, pol::Vector, sense="inf")
    if pol != nothing
        f = [p*one(Polynomial{true,Float64}) for p in pol]
        MOM.set_objective(M,idx,f, sense)
    end
end

#----------------------------------------------------------------------
function objective_ncl(M)
    B = M[:basis]
    @objective(M.model, Min,  sum(sum(M[:moments][M[:index][B[i]^2],k] for i in 1:length(B)) for k in 1:M[:nu]))
end
#----------------------------------------------------------------------
function objective_tv(M)
    objective(M,collect(1:M[:nu]), ones(M[:nu]),"inf")
end
