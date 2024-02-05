export polar_ideal, polar_minimize, preordering

#----------------------------------------------------------------------
"""
```julia
j = polar_ideal(f, [h1, h2, ...], [g1, g2, ...], X)
```
Compute generators of the polar ideal associated with `f, [h1, h2, ...], [g1, g2, ...]` (equality constraints hi == 0 and the sign constraints gi >= 0) with respect to the variable vector X.

f, hi, gi should be polynomials in variables X.
"""
function polar_ideal(f, h::Vector, g::Vector, X)
    n = length(X)
    r = length(h)

    if r+1 > n
        println("Too many equality constraints: polar ideal gives nothing.")
        return nothing
    end

    # Converting everithing in Polynomial type
    f = one(f)*f
    h = one(f)*h
    g = one(f)*g

    # First work without nonneg constraints (to avoid emptysequence problems)
    partial_matrix , eq = zero_iteration(X, f, h)

    #Adding the nonneg constraints
    eq = nonneg_iterations(X, partial_matrix, eq, g, r)

    # Return equations without zeros and without repetitons
    return setdiff([h...,eq...], 0)
end

#----------------------------------------------------------------------
function zero_iteration(X, f, h)
    n = length(X)
    r = length(h)

    # Condider the emptysequence case separately because  differentiate([f; []], X) doesn't work
    if h == []
        partial_matrix = differentiate([f],X)
        eq = transpose(partial_matrix)
    else
        partial_matrix = differentiate([f; h], X)
        eq = Polynomial[]
        for j in combinations(1:n, r+1)
            push!(eq, det(partial_matrix[:,j]))
        end
    end
    return partial_matrix, setdiff(eq, 0)
end

#----------------------------------------------------------------------
function nonneg_iterations(X, partial_matrix, eq, g, r)
    n = length(X)
    for k in 1:length(g)
        for gk in collect(combinations(g,k))
            gk = reshape(gk,length(gk),1)
            # if there are too many rows the det equations are trivial
            if k < n-r
                matrix = [partial_matrix ; differentiate(gk,X)]
                for j in collect(combinations(1:n, r+k+1))
                    temp = det(matrix[:,j])
                    #if det is a number != 0 then ideal is (1); if det = 0 then skip.
                    if temp isa Number
                        if temp == 0
                            continue
                        else
                            gk = [1]
                            break
                        end
                    else
                        gk = vcat(gk, temp)
                    end
                end
            end
            eq = kron(eq, gk)
        end
    end
    return setdiff(eq, 0)
end

#----------------------------------------------------------------------
"""
```julia
v, M = polar_minimize(f, [h1, h2, ...], [g1, g2, ...], X, degree_approx, optimizer)
```
Compute the infimum of the f (equality constraints hi == 0 and the sign constraints gi >= 0) on the moment side.
It does that calling minimize(...), replacing the equality constraints  [h1, h2, ...] with generators of the polar ideal.

f, hi, gi should be polynomials in the variables X.

If the problem is feasible and has minimizers, it outputs
  - v: the minimum value
  - M: the moment model of type MOM.Model
"""
function polar_minimize(f, h::Vector, g::Vector, X, d::Int64)
    j = polar_ideal(f, h, g, X)
    return optimize(:inf, f, j, g, X, d)
end

#----------------------------------------------------------------------
"""
```julia
pg = preordering([g1, g2, ...])
```
Compute all the product, without repetitions, of the gi's i.e. computes the generators (as a quadratic module) of the preordering O(g1, g2, ...).

"""
function preordering(g::Vector)
    pg = [1];
    for p in collect(powerset(g,1,length(g)))
        reshape(p,length(p),1)
        pg = vcat(pg, prod(p))
    end
    return setdiff(pg, 1)
end
