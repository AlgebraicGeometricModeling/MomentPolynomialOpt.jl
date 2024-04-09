export MOM, moments
#import JuMP: optimize!, objective_value


"""
```
moments(M, X, d::Int, symb::Symbol...)
```

Define variables for the moments of the monomials of degree less than d in the variables  X.
If
 - `symb == PSD`, non-negativity constraints are added;
 - `symb == PRB`, unit mass and non-negativity constraints are added.

"""
function moments(M, X, d::Int, symb::Symbol...)
    mu = MOM.moment_variables(M, X, d)
    for arg in symb
        if arg == :PRB
            MOM.add_constraint_nneg(M,mu)
            MOM.add_constraint_unitmass(M, mu)
        elseif arg == :PSD
            MOM.add_constraint_nneg(M,mu)
        end
    end
    return mu 
end


module MOM


using DynamicPolynomials, JuMP, Dualization, LinearAlgebra

import MomentPolynomialOpt: MPO, Moments, moments, mmt

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

function moment_variables(M, B)
    v = @variable(M, [1:length(B)])
    mu = Moments(v,B)
    if haskey(M.obj_dict, :mu)
        push!(M[:mu],mu)
    else
        M[:mu] = [mu]
    end
    return mu
end

function moment_variables(M, X, d::Int)
    moment_variables(M, monomials(X,0:d))
end

#----------------------------------------------------------------------
function Model(optimizer = MPO[:optimizer])

    if optimizer == nothing
        @error "No optimizer set; see mpo_optimizer"
    else
        M = JuMP.Model(Dualization.dual_optimizer(optimizer))
    end
    
    M[:type] = :moment

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

    M = MOM.Model(optimizer)

    mu = moments(M, X, 2*d, :PRB)
    
    for h in H MOM.add_constraint_zero(M, mu, h) end
    for g in G MOM.add_constraint_nneg(M, mu, g) end
    
    @objective(M, Min, mmt(mu,f))

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

    M = MOM.Model(optimizer)
    
    mu = moments(M, X, 2*d, :PRB)

    wobj = false
    for c in C
        if c[2] == "inf" || c[2] == "min"
            @objective(M, Min, mmt(mu,c[1]))
            wobj = true
        elseif c[2] == "sup" || c[2] == "max"
            MOM.set_objective(M, "sup", mu, c[1])
            wobj = true
        elseif c[2] == "=0"
            MOM.add_constraint_zero(M, mu, c[1])
        elseif c[2] == ">=0"
            MOM.add_constraint_nneg(M, mu, c[1])
        elseif c[2] == "<=0"
            MOM.add_constraint_nneg(M, mu, -c[1])
        end
    end
    if !wobj
        @objective(M, Max, mmt(mu, 1)) #MOM.set_objective(M, "sup", one(C[1][1]), mu)
    end

    return M
    
end

end  #module MOM

