export get_series, get_optimizers, get_measure

using MultivariateSeries

#=
function dualize(p)
    MultivariateSeries.dual(p)
end
=#

#----------------------------------------------------------------------
"""
```
get_series(M)
```
Return the vector of ``\\nu``=`M[:nu]` series of optimal moments of the optimized moment program `M`.
"""
function get_series(M::JuMP.Model)

    if haskey(M.obj_dict,:type)
        if M[:type] == :polynomial
            s = get_series_dual.(M[:mu])
            return (length(s) == 1 ? s[1] : s) 
        elseif M[:type] == :moment
            s = get_series_primal.(M[:mu])
            return (length(s) == 1 ? s[1] : s) 
        end
    else
        if haskey(M.obj_dict,:constraints) && haskey(M.obj_dict,:monomials) 
            return [MultivariateSeries.dual(M[:monomials]'*JuMP.dual.(M[:constraints]))]
        else
            n = length(M[:monomials])
            cstr = JuMP.all_constraints(M[:dual], AffExpr, MOI.EqualTo{Float64})
            s  = [MultivariateSeries.series([M[:monomials][i]=>-JuMP.dual(cstr[n*(k-1)+i])
                                             for i in 1:n])
                  for k in 1:M[:nu]]
            
            return s
            
            [MultivariateSeries.series([M[:monomials][i]=>JuMP.value(M[:moments][k,i])
                                        for i in 1:length(M[:monomials])])
             for k in 1:M[:nu]]
            
        end
    end
end
    
function get_series_primal(s :: Moments)
        MultivariateSeries.series([s.basis[i] => JuMP.value(s.values[i]) for i in 1:length(s.basis)])
end
    
function get_series_dual(s :: Moments)
        MultivariateSeries.series([s.basis[i] => JuMP.dual(s.values[i]) for i in 1:length(s.basis)])
end

#----------------------------------------------------------------------
"""
```
get_optimizers(M, , t::Int64 = Inf)
```
Return the optimal points  of the moment program `M`, using moments of degree <=t
(default: twice the order of the relaxation minus 2)

```julia
get_optimizers(M)

[1.41421 1.73205; 1.41421 1.41421; 1.41421 -1.73205]
```
"""
function get_optimizers(M::JuMP.Model)
    s = get_series(M)
    t = maxdegree(s)-1
    w, Xi = MultivariateSeries.decompose(truncate(s, t));
    Xi
end

#----------------------------------------------------------------------
"""
```
w, Xi = get_measure(M)
```
Return the approximation of the moment sequence ``\\mu``
truncated to moments of degree <= t (default: twice the order of the relaxation minus 2),
as weighted sum of Dirac measures: ``\\sum_{k=1}^{r} \\omega_k \\delta_{\\xi_k}`` where

- `w` is the vector of weights of the Dirac measures.
- `Xi` is matrix of ``n\\times r`` support points of the corresponding Dirac measures. The column `Xi[:,i]` is the support point ``\\xi_{i}`` of the ith Dirac measure and its weights is `w[i]`.

```julia
w, Xi = get_measure(M)

([0.1541368146508854, 0.5889741915171074, 0.256888993597116], [1.4142135624216647 1.414213562080608 1.4142135620270329; -1.732052464639053 1.4141771454788292 1.7319839273833693])
```
"""
function get_measure(M::JuMP.Model)

    s = get_series(M)
    t = maxdegree(s)-1
    w, Pts = MultivariateSeries.decompose(truncate(s,t));

    return w, Pts
end

function get_measure(M::JuMP.Model,  e::Float64)

    s = get_series(M)
    t = maxdegree(s)-1
    w, Pts = MultivariateSeries.decompose(truncate(s,t), MultivariateSeries.eps_rkf(e));
    return w, Pts
end





