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
