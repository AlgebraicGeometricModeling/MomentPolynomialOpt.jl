export MinEllipsoidModel

lower_triangular(P,s=1) = [P[i, j] for i = s:size(P, 1) for j = s:i]

# Create a Sos + LogDetCone Model for the minimum ellipsoid
function MinEllipsoidModel(pt, H, G, X, d)

    M = JuMP.Model()

    Mh = [monomials(X,0:2*d-maxdegree(h)) for h in H]

    Lg = [monomials(X,0:d-Int64(ceil(maxdegree(g)/2))) for g in G]
    L0 = monomials(X,0:d)


    a = [@variable(M, [1:length(Mh[i])], base_name="a$i") for i in 1:length(H)]

    Q = [@variable(M, [1:length(Lg[i]), 1:length(Lg[i])], PSD, base_name="Q$i") for i in 1:length(G)]
    Q0 = @variable(M, [1:length(L0), 1:length(L0)], PSD, base_name="Q0") 

    n = length(pt)
    
    @variable(M, Z[1:n+1,1:n+1], PSD)

    @variable(M, t)
    
    M[:a] = a
    M[:Q] = Q
    M[:Q0]= Q0

    # Constraints
    @constraint(M, [t; 1; lower_triangular(Z,2)] in MOI.LogDetConeTriangle(n))

    p = [-1. ; pt ]
    
    P = 1 - p'*Z*p -
        sum( H[i]*sum(a[i][j]*Mh[i][j] for j in 1:length(Mh[i])) for i in 1:length(H)) -
        sum( G[i]*(Lg[i]'*Q[i]*Lg[i]) for i in 1:length(G) ) -
        L0'*Q0*L0

    for c in coefficients(P)
        @constraint(M, c == 0)
    end

    @objective(M, Max, t)
    return M
    
end

export min_ellipsoid


"""
Compute the minimal ellipsoid enclosing the point P of the semi-algebraic set defined by H[j]=0, G[k]>=0

It outputs the center c and the matrix U which columns are the principal axes of the ellipsoid.
"""
function min_ellipsoid(P, H, G, X, d, optimizer=MMT[:optimizer])

    M = MinEllipsoidModel(P,H,G,X,d)
    set_optimizer(M,optimizer)
    optimize!(M)
    Z = M[:Z]; 

    D = value.(Z[2:end,2:end]); 
    c = inv(D)*value.(Z[1,2:end]);
    
    U, S, V = svd(Matrix(D))

    for i in 1:size(U,2) U[:,i]/= sqrt(S[i]) end
    
    return c, U, M
    
end

function min_ellipsoid(x, optimizer=MMT[:optimizer])

    M = JuMP.Model()

    n = size(x,1)
    #@variable(model, z[1:n])
    @variable(M, Z[1:n+1,1:n+1], PSD)
    @variable(M, t)
    #@variable(M,s)
    
    #@constraint(M, Z in PSDCone())
    #@constraint(M, [s z'; z Z] in PSDCone())
    @constraint(M, [t; 1; lower_triangular(Z,2)] in MOI.LogDetConeTriangle(n))
    
    @inbounds for i=1:size(x,2)
        p = [-1., x[:,i]...]
        @constraint(M, p'*Z*p <= 1)
    end

    @objective(M, Max, t)

    set_optimizer(M,optimizer)
    # solve and query solution
    optimize!(M)
    # println(termination_status(M))
    # println(objective_value(M ))
    
    D = value.(Z[2:end,2:end]); 
    c = inv(value.(Z[2:end,2:end]))*value.(Z[1,2:end]);

    U, S, V = svd(Matrix(D))

    for i in 1:size(U,2) U[:,i]/= sqrt(S[i]) end
    
    return c, U, M
    
end
