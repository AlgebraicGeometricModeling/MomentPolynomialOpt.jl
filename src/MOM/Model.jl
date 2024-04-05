export MOM
import JuMP: optimize!, objective_value

module MOM


using DynamicPolynomials, JuMP, Dualization

# import MathOptInterface
# const MOI = MathOptInterface

using LinearAlgebra

import MomentPolynomialOpt: MPO, Moments, mmt

export Model

#----------------------------------------------------------------------
convert_Float64 = function(pol)
    if typeof(pol) != Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Int64}
        return dot(Float64.(coefficients(pol)),monomials(pol))
    else
        return pol
    end
end
#----------------------------------------------------------------------

include("constraints.jl")
include("objective.jl")

function moment_variables(M, B, cstr = :PSD)
    v = @variable(M, [1:length(B)])
    mu = Moments(v,B)
    if haskey(M.obj_dict, :mu)
        push!(M[:mu],mu)
    else
        M[:mu] = [mu]
    end
    if cstr == :PSD
        add_constraint_nneg(M, mu)
    end
    return mu
end

function moment_variables(M, X, d::Int, cstr = :PSD)
    moment_variables(M, monomials(X,0:d), cstr)
end

#----------------------------------------------------------------------
function Model(X, d, optimizer = MPO[:optimizer])

    if optimizer == nothing
        @error "No optimizer set; see mpo_optimizer"
    else
        M = JuMP.Model(Dualization.dual_optimizer(optimizer))
    end
    
    M[:type] = :moment
    M[:variables] = X
    M[:degree] = d
    
    mu = MOM.moment_variables(M, X, 2*d, :PSD) 
    # MOM.add_constraint_nneg(M, mu)

    return M
end

"""
```
M = MOM.Model( `sense`, f, H, G, X, d)
```
Construct the Moment Program in the variables X of order d.
   - `sense` == "inf" or "sup"
   - `f` polynomial objective function
   - `H =[h1, h2, ...]` array of polynomial equality constraints (can be empty)
   - `G =[g1, g2, ...]` array of non-negativity constraints (can be empty)
   - `X` is the vector of variables
   - `d` is the order of the moment relaxation.
"""
function Model(sense::Symbol, f, H, G, X, d, optimizer = MMT[:optimizer])

    M = MOM.Model(X,d,optimizer)

    mu = M[:mu][1]

    MOM.add_constraint_unitmass(M, mu)
    
    for h in H MOM.add_constraint_zero(M, h, mu) end
    for g in G MOM.add_constraint_nneg(M, g, mu) end
    
    MOM.set_objective(M, "inf", f, mu)

    return M
    
end



"""
```
M = MOM.Model(C, X, d)
```
Construct the Moment Program where
   - C is a vector of pairs (f, sense ) of objective or constraints where f is a polynomial and sense is "inf", "min", "sup", "max", ">=0", "<=0", "=0", or an interval 
   - `X` is the vector of variables
   - `d` is the order of the moment relaxation.
"""
function  Model(C::Vector, X, d::Int64, optimizer = MMT["optimizer"]; kwargs...)

    M = MOM.Model(X,d,optimizer)
    
    mu = M[:mu][1]   

    MOM.add_constraint_unitmass(M, mu)

    wobj = false
    for c in C
        if c[2] == "inf" || c[2] == "min"
            MOM.set_objective(M, "inf", c[1], mu)
            wobj = true
        elseif c[2] == "sup" || c[2] == "max"
            MOM.set_objective(M, "sup", c[1], mu)
            wobj = true
        elseif c[2] == "=0"
            MOM.add_constraint_zero(M, c[1], mu)
        elseif c[2] == ">=0"
            MOM.add_constraint_nneg(M, c[1], mu)
        elseif c[2] == "<=0"
            MOM.add_constraint_nneg(M,-c[1], mu)
        end
    end
    if !wobj
        MOM.set_objective(M, "sup", one(C[1][1]), mu)
    end

    return M
    
end

end  #module MOM

