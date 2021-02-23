export annihilator, is_aprox_zero, is_zero;

#------------------------------------------------------------------------
is_aprox_zero(eps) = function(d)
  abs(d) < eps
end

# function is_zero(d::T) where {T}
#   d == zero(T)
# end

function is_zero(d::Float64)
    is_aprox_zero(0.000001)(d)
end
#----------------------------------------------------------------------


function orthoproj(sigma,f,P)
   r=f
    for i in 1:length(P)
        n2 = LinearAlgebra.dot(sigma,P[i]*P[i])
        r -= LinearAlgebra.dot(sigma,r*P[i])/n2*P[i]
   end
   r
end


#----------------------------------------------------------------------
"""
```
K,I,P,B = annihilator(sigma::Series{R}, L, eqzero = is_zero)
```
Compute a graded basis of the annihilator of the positive series ``σ``
in the space spanned by the list of monomials L.

It outputs:

  * K the graded basis
  * I the initial monomials of K
  * P the orthogonal basis
  * B the monomial basis

"""
function annihilator(sigma, L, eps = 1.e-4)
    X = variables(sigma)
    dA = deg(L[length(L)])

    P = []
    B = []
    K = []
    I = []

    for m in L
        p  = orthoproj(sigma,m,P)

        n2 = dot(sigma, p*p)

        if abs(n2) < eps
            #println("n2: ",m," ", n2, " ", eps)
            nw = true
            for lt in I
                if divides(lt,m)
                    nw = false
                    break
                end
            end
            if nw
                push!(K,p)
                push!(I,m)
            end
        else
            push!(P,p)
            push!(B,m)
        end
    end

    K,I,P,B
end
#------------------------------------------------------------------------
