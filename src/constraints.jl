export constraint_zero, constraint_nneg, constraint_unitmass, constraint_moments

#----------------------------------------------------------------------
function constraint_zero(M::MOM.Model, idx::Int64, eqs...)
    for e in eqs
        MOM.add_constraint_zero(M,[idx],[e])
    end
end

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
function constraint_nneg(M::MOM.Model, idx::Int64, eqs...)
    for e in eqs
        MOM.add_constraint_nneg(M,[idx],[e])
    end
end

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

function constraint_moments(M::MOM.Model, moments::Vector, idx::Int64, p)
    for c in moments
        MOM.add_constraint_moment(M, c[2], [idx], [p*c[1]*one(Polynomial{true,Float64})])
    end
end

function constraint_moments(M::MOM.Model, moments::Vector, idx::Vector{Int64}, P::Vector)
    for c in moments
        MOM.add_constraint_moment(M, c[2], idx, P*c[1]*one(Polynomial{true,Float64}))
    end
end

#----------------------------------------------------------------------
function constraint_unitmass(M::MOM.Model)
    MOM.add_constraint_moment(M, 1, one(Polynomial{true,Float64}))
end

function constraint_unitmass(M::MOM.Model, idx::Int64, p)
    MOM.add_constraint_moment(M, 1, [idx], [p])
end

function constraint_unitmass(M::MOM.Model, idx::Vector{Int64}, P::Vector)
    MOM.add_constraint_moment(M, 1, idx, P)
end

#----------------------------------------------------------------------
