export MaxEigenModel

function MaxEigenModel(f, H::AbstractVector, G::AbstractVector, X, d::Int64, optimizer = MMT[:optimizer])

    M = JuMP.Model(optimizer)
    
    M[:type] = :polynomial
    
    MP = [monomials(X,0:2*d-maxdegree(h)) for h in H]
    M[:mP] = MP
    if length(G)>0
        MQ = [monomials(X,0:d-Int64(ceil(maxdegree(g)/2))) for g in G]
    else
        MQ=[]
    end
    
    M[:mQ] = MQ
    L0 = monomials(X,0:d)
    M[:monomials] = L0

    @variable(M, lambda)

    P = [@variable(M, [1:length(MP[i])], base_name="P$i") for i in 1:length(H)]

    if length(G) > 0
        Q = [@variable(M, [1:length(MQ[i]), 1:length(MQ[i])], PSD, base_name="Q$i") for i in 1:length(G)]
    else
        Q = []
    end
    
    n = length(L0)
    Q0 = @variable(M, [1:n, 1:n], Symmetric, base_name="Q0") 

    # Positivity constraint: Q - lambda*I >=0
    @constraint(M, Q0 - diagm(fill(lambda,n)) in PSDCone())
        
    M[:P] = P
    M[:Q] = Q
    M[:Q0]= Q0

    R = f - L0'*Q0*L0

    for i in 1:length(H)
        R -= H[i]*sum(P[i][j]*MP[i][j] for j in 1:length(MP[i]))
    end
    
    for i in 1:length(G)
        R -= G[i]*(MQ[i]'*Q[i]*MQ[i])
    end 
   
#    for c in coefficients(R)
#        @constraint(M, c == 0)
#    end

    M[:mu] = [ Moments([@constraint(M, c == 0) for c in coefficients(R)], monomials(R)) ]
    
    @objective(M, Max, lambda)

    return M
end
    
