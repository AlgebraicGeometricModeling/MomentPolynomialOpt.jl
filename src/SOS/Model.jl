
#export SOSModel
export SOS

module SOS

import MomentTools: MMT


using JuMP, DynamicPolynomials

export Model

function Model(sense::Symbol, f, H, G, X, d, optimizer = MMT[:optimizer])

    if optimizer == nothing
        M = JuMP.Model()
    else
        M = JuMP.Model(with_optimizer(optimizer))
    end
    
    M[:type] = :polynomial
    
    Mh = [monomials(X,0:2*d-maxdegree(h)) for h in H]

    Lg = [monomials(X,0:d-Int64(ceil(maxdegree(g)/2))) for g in G]
    L0 = monomials(X,0:d)

    @variable(M, lambda)

    a = [@variable(M, [1:length(Mh[i])], base_name="a$i") for i in 1:length(H)]

    Q = [@variable(M, [1:length(Lg[i]), 1:length(Lg[i])], PSD, base_name="Q$i") for i in 1:length(G)]
    Q0 = @variable(M, [1:length(L0), 1:length(L0)], PSD, base_name="Q0") 
        
    M[:a] = a
    M[:Q] = Q
    M[:Q0]= Q0

    if in(sense, [:Min, :min, :Inf, :inf])
        
        P = f - lambda - L0'*Q0*L0
        for i in 1:length(H) P -= H[i]*sum(a[i][j]*Mh[i][j] for j in 1:length(Mh[i])) end
        for i in 1:length(G) P -= G[i]*(Lg[i]'*Q[i]*Lg[i]) end

        @objective(M, Max, lambda)

    else
        
        P = lambda - f - L0'*Q0*L0
        for i in 1:length(H) P -= H[i]*(a[i]'*Mh[i]) end
        for i in 1:length(G) P -= G[i]*(Lg[i]'*Q[i]*Lg[i]) end
           

        @objective(M, Min, lambda)
    end

    for c in coefficients(P)
        @constraint(M, c == 0)
    end


    return M
end

#=
function optimize(sense::Symbol, f, H, G, X, d, optimizer = MMT[:optimizer])
    M = SOS.Model(sense,f,H,G,X,d, optimizer)
    optimize!(M)

    return objective_value(M), M
end

function minimize(f, H, G, X, d , optimizer = MMT[:optimizer])
   return SOS.optimize(:inf, f,H,G,X,d, optimizer)
end

function maximize(f, H, G, X, d , optimizer = MMT[:optimizer])
   return SOS.optimize(:sup, f,H,G,X,d, optimizer)
end
=#
end #module SOS


