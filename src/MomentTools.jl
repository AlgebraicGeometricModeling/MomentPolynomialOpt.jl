module MomentTools

using JuMP
using DynamicPolynomials
using MultivariateSeries

export Seq, seq

mutable struct Seq{T} val::T end

function seq(args...)
    if length(args)>1
        Seq([args...])
    else
        Seq(args[1])
    end
end

function DynamicPolynomials.MonomialVector(V::Vector{PolyVar{true}}, rg::Seq)
     L = DynamicPolynomials.Monomial{true}[]
     for i in rg.val
         append!(L, DynamicPolynomials.monomials(V,i))
     end
     L
end


include("mommodel.jl")
include("constraints.jl")
include("optimize.jl")


export MOM, getseries, getminimizers, getmeasure

end
