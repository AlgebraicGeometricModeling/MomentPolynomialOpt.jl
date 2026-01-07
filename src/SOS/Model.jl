export SOS

module SOS

import MomentPolynomialOpt: MPO, Moments

using JuMP, DynamicPolynomials

export Model

function Model end

function Model(sense::Symbol, f, H, G, X, d, optimizer = MPO[:optimizer])

    if optimizer == nothing
        M = JuMP.Model()
    else
        M = JuMP.Model(optimizer)
    end
    
    M[:type] = :polynomial
    
    Lh = [monomials(X,0:2*d-maxdegree(h)) for h in H]
    Lg = [monomials(X,0:d-Int64(ceil(maxdegree(g)/2))) for g in G]
    L0 = monomials(X,0:d)

    @variable(M, lambda)

    p = [@variable(M, [1:length(Lh[i])], base_name="p$i") for i in 1:length(H)]

    Q =[@variable(M, [1:length(L0), 1:length(L0)], PSD, base_name="Q0")]
    append!(Q,[@variable(M, [1:length(Lg[i]), 1:length(Lg[i])], PSD, base_name="Q$i") for i in 1:length(G)])
    
        
    M[:p] = p
    M[:Q] = Q

    if in(sense, [:Min, :min, :Inf, :inf])
        
        P = f - lambda - L0'*Q[1]*L0
        for i in 1:length(H)
            P -= H[i]*sum(p[i][j]*Lh[i][j] for j in 1:length(Lh[i]))
        end
        for i in 1:length(G)
            P -= G[i]*(Lg[i]'*Q[i+1]*Lg[i])
        end

        @objective(M, Max, lambda)

    else
        
        P = lambda - f - L0'*Q[1]*L0
        for i in 1:length(H) P -= H[i]*(p[i]'*Lh[i]) end
        for i in 1:length(G) P -= G[i]*(Lg[i]'*Q[i+1]*Lg[i]) end
           

        @objective(M, Min, lambda)
    end

    M[:mu] = [ Moments([@constraint(M, c == 0) for c in coefficients(P)], monomials(P)) ]

    return M
end

#=
function optimize(sense::Symbol, f, H, G, X, d, optimizer = MPO[:optimizer])
    M = SOS.Model(sense,f,H,G,X,d, optimizer)
    optimize!(M)

    return objective_value(M), M
end

function minimize(f, H, G, X, d , optimizer = MPO[:optimizer])
   return SOS.optimize(:inf, f,H,G,X,d, optimizer)
end

function maximize(f, H, G, X, d , optimizer = MPO[:optimizer])
   return SOS.optimize(:sup, f,H,G,X,d, optimizer)
end
=#
end #module SOS


