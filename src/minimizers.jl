export get_series, get_minimizers, get_measure, dualize

function dualize(p)
    MultivariateSeries.dual(p)
end

#----------------------------------------------------------------------
"""
```
get_series(M)
```
Return the vector of ``\\nu``=`M[:nu]` series of optimal moments of the optimized moment program `M`.
"""
function get_series(M::JuMP.Model)

    if haskey(M.obj_dict,:constraints) && haskey(M.obj_dict,:monomials) 
        return [MultivariateSeries.dual(M[:monomials]'*JuMP.dual.(M[:constraints]))]
    else
        n = length(M[:monomials])
        cstr = JuMP.all_constraints(M[:dual], AffExpr, MOI.EqualTo{Float64})
        s  = [series([M[:monomials][i]=>-JuMP.dual(cstr[n*(k-1)+i])
                  for i in 1:n])
              for k in 1:M[:nu]]

        return s
    
        [series([M[:monomials][i]=>JuMP.value(M[:moments][k,i])
                 for i in 1:length(M[:monomials])])
         for k in 1:M[:nu]]
        
    end
end

#----------------------------------------------------------------------
"""
```
get_minimizers(M, , t::Int64 = 2*M[:degree]-1)
```
Return the minimizer points  of the optimized moment program `M`, using moments of degree <=t
(default: twice the order of the relaxation minus 2)

```julia
get_minimizer(M)

[1.41421 1.73205; 1.41421 1.41421; 1.41421 -1.73205]
```
"""
function get_minimizers(M::JuMP.Model, t::Int64 = 2*M[:degree]-1)
    s = get_series(M)
    w, Xi = MultivariateSeries.decompose(truncate(s[1], t));
    Xi
end

#----------------------------------------------------------------------
"""
```
w, Xi = get_measure(M, t::Int64 = 2*M[:degree]-1 ,lambda = [(-1)^(k-1) for k in 1:M[:nu]])
```
Return the approximation of the moment sequence ``\\sum_{i=1}^{\\nu} \\lambda_i \\mu_i``
truncated to moments of degree <= t (default: twice the order of the relaxation minus 2),
as weighted sum of Dirac measures: ``\\sum_{k=1}^{r} \\omega_k \\delta_{\\xi_k}`` where

- `w` is the vector of weights of the Dirac measures.
- `Xi` is matrix of ``n\\times r`` support points of the corresponding Dirac measures. The column `Xi[:,i]` is the support point ``\\xi_{i}`` of the ith Dirac measure and its weights is `w[i]`.

```julia
w, Xi = get_measure(M)

([0.1541368146508854, 0.5889741915171074, 0.256888993597116], [1.4142135624216647 1.414213562080608 1.4142135620270329; -1.732052464639053 1.4141771454788292 1.7319839273833693])
```
"""
function get_measure(M::JuMP.Model,
                     t::Int64 = 2*M[:degree]-1,
                     lambda::Vector = [(-1)^(k-1) for k in 1:get(M.obj_dict,:nu,1)])

    s = get_series(M)
    w, Pts = MultivariateSeries.decompose(truncate(s[1], t));
    nu = get(M.obj_dict,:nu,1)
    if nu>1
        for k in 2:nu
            c, Xi = MultivariateSeries.decompose(truncate(s[k], t))
            w = vcat(w, c*lambda[k])
            Pts= hcat(Pts,Xi)
        end
    end

    return w, Pts
end

function get_measure(M::JuMP.Model, e::Float64,
                     t::Int64 = 2*M[:degree]-2,
                     lambda::Vector = [(-1)^(k-1) for k in 1:get(M.obj_dict,:nu,1)])

    s = get_series(M)
    w, Pts = MultivariateSeries.decompose(truncate(s[1],t), MultivariateSeries.eps_rkf(e));
    nu = get(M.obj_dict,:nu,1)
    if nu>1
        for k in 2:nu
            c, Xi = MultivariateSeries.decompose(truncate(s[k],t), MultivariateSeries.eps_rkf(e))
            w = vcat(w, c*lambda[k])
            Pts= hcat(Pts,Xi)
        end
    end
    return w, Pts
end





