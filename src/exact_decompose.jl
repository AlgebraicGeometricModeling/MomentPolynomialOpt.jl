using MultivariateSeries

export exact_decompose


function Q_matrix(g, n)
    Q  = fill(zero(coefficients(g)[1]),n,n)
    for t in g
        i = maxdegree(t)
        c = coefficient(t)
        if i <= n-1
            if i == 0
                Q[1,1]= c
            else
                Q[1,i+1] = c/2
            end
        else
            if i == 2n-2
                Q[n,n]= c
            else
                Q[i+2-n,n] = c/2
            end
        end
    end
    return Symmetric(Q)
end

function H_matrices(n, d)
    [
        [(i+j==k+2 ? 1//1 : 0//1 ) for i in 1:n, j in 1:n ]
        for k in 0:d
    ]
end

function H_proj(Q, H)
    sum(h*tr(h*Q)/tr(h*h) for h in H) 
end


"""
Decompose g>0 at the real roots of f as a Sum of Squares and a multiple of f with rational 
coefficients.

  - Q rational matrix which is positive definite if g>0 at the real roots of f
  - q rational factor 

such that g = x^t Q x + q*f.

The options are:

   - rounding = (Int64)k upper bound on the precision used in the rounding step
   - optimizer = optimizer tool used for solving the SDP problem
   - verbose = true or false to specify if information during the computation is printed or not

"""
function exact_decompose(g, f; rounding = 10000, optimizer = MMT[:optimizer], verbose = false)
    n = maxdegree(f)
    dg = maxdegree(g)
    dq = n-dg

    X  = variables(f)
    L  = monomials(X,0:2n-2)
    B  = monomials(X,0:n-1)
    Lq = monomials(X,0:dq)

    cf = MultivariateSeries.matrixof([f],L)
    cg = MultivariateSeries.matrixof([g],L)

    M = JuMP.Model()

    
    set_optimizer(M, optimizer)

    # The variables are lambda, the symmetric matrix Q and the coefficients of q
    @variable(M, lambda)
    # @variable(M, Q[1:n, 1:n], PSD)
    @variable(M, Q[1:n, 1:n], Symmetric)
    @variable(M, q[1:dq+1])


    # Positivity constraint: Q - lambda*I >=0
    @constraint(M, Q - lambda*Matrix(I,n,n) in PSDCone())

    # The constraints g = x^t Q x + q*f
    for k in 1:length(L)
        s = min(n,k)
        @constraint(M,
                    sum(Q[s+1-i,k+i-s] for i in 1:min(k,2n-k))
                    + sum(q[i]*cf[k+1-i] for i in 1:min(dq,k))
                    - cg[1,k] == 0)
    end    

    # Objectiv function = lambda to maximize the smallest eigenvalmue of Q

    @objective(M, Max, lambda)

    JuMP.optimize!(M)

    lambda = JuMP.value.(lambda)

    if verbose 
        println("Minimal eigenvalue: ", lambda)
    end
    
    Qopt = JuMP.value.(Q)
    qopt = JuMP.value.(q)
    gopt = B'*Qopt*B + (qopt'*Lq)*f

    H = H_matrices(n,2n-2)

    popt = - (qopt'*Lq)*f + g
    Qpopt = Q_matrix(popt,n)

    epol = g - (qopt'*Lq)*f - B'*Qopt*B

    rho = norm(epol)
    
    # Ev = eigen(Qopt).values
    # println("Eigen: ", Ev)
    
    # Rounding
    #delta=0.99*(min(Ev...)-rho)/(n+ sqrt(n)*(n-1)*norm(f))
    delta=0.99*(lambda-rho)/(n+ sqrt(n)*(n-1)*norm(f))
    #println("delta: ", delta)

    k = Int64(ceil(log10(1/delta)))
    k = min(k, rounding)

    if verbose 
        println("Rounding precision: $k digit(s)")
    end

    Qr = Int64.(round.(Qopt*10^k; digits = 0))//10^k
    qr = Int64.(round.(qopt*10^k; digits = 0))//10^k
 
    
    
    # Projection
    qexact = (qr'*Lq)
    p  = -qexact*f + g
    Qp = Q_matrix(p,n)
    Qexact = Qr + H_proj(Qp-Qr, H)
    
    
    if verbose 
        # Certificate
        delta = B'*Qexact*B+qexact*f-g
        # println("delta: ",delta)
        @assert delta == 0
    end
    
    Qexact, qexact
end



