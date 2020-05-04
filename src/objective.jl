export objective, objective_ncl, objective_tv

"""
 Add as objective function the linear functional associated to the polynomial pol to minimize.
"""
function objective(M::MOM.Model, pol, sense="inf")
    if pol != nothing
        MOM.set_objective(M,pol,sense)
    end
end

function objective(M::MOM.Model, pol::Vector, sense="inf")
    if pol != nothing
        MOM.set_objective(M, collect(1:M[:nu]), pol, sense)
    end
end

function objective(M::MOM.Model, idx:: Int64, pol, sense="inf")
    if pol != nothing
        MOM.set_objective(M, [idx], [pol], sense)
    end
end


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
    @objective(M.model, Min,  sum(sum(M[:moments][k,M[:index][B[i]^2]] for i in 1:length(B)) for k in 1:M[:nu]))
end
#----------------------------------------------------------------------
function objective_tv(M)
    objective(M,collect(1:M[:nu]), ones(M[:nu]),"inf")
end

