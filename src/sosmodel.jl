export SOS

module SOS

using JuMP, DynamicPolynomials

import MomentTools: MMT

function Model(sense::Symbol, f, H, G, X, d)

    M = JuMP.Model()
    
    if haskey(MMT,:optimizer,)
        set_optimizer(M, MMT[:optimizer])
    end


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

    if sense == :Min
        
        P = f - lambda -
            sum( H[i]*sum(a[i][j]*Mh[i][j] for j in 1:length(Mh[i])) for i in 1:length(H)) -
            sum( G[i]*(Lg[i]'*Q[i]*Lg[i]) for i in 1:length(G) ) -
            L0'*Q0*L0

        @objective(M, Max, lambda)

    else
        
        P = lambda - f -
            sum( H[i]*sum(a[i][j]*Mh[i][j] for j in 1:length(Mh[i])) for i in 1:length(H)) -
            sum( G[i]*(Lg[i]'*Q[i]*Lg[i]) for i in 1:length(G) ) -
            L0'*Q0*L0

        @objective(M, Min, lambda)
    end

    for c in coefficients(P)
        @constraint(M, c == 0)
    end


    return M
end


function optimize(sense::Symbol, f, H, G, X, d)
    M = SOS.Model(sense,f,H,G,X,d)
    optimize!(M)

    return objective_value(M), M
end

end #mod


