export Moments, mmt

mutable struct Moments{M}
    values::AbstractVector
    basis::AbstractVector{M}
    index::Dict{M,Int}
end

function Moments(v::Vector, B)
    idx= Dict{typeof(B[1]),Int}()
    for i in 1:length(B)
        idx[B[i]]=i
    end
    Moments(v,B,idx)
end

function mmt(mu::Moments, c::Real)
    mu.values[1]*c
end

function mmt(mu::Moments, x::Variable)
    mmt(mu, monomial(x))
end

function mmt(mu::Moments, m::Monomial)
    mu.values[mu.index[m]]
end

function mmt(mu::Moments, p::AbstractPolynomial)
    sum(t.coefficient*mu.values[mu.index[t.monomial]] for t in terms(p))
end

function moment_variables(M, B, cstr = :PSD)
    v = @variable(M, [1:length(B)])
    mu = Moments(v,B)
    if haskey(M.obj_dict, :mu)
        push!(M[:mu], mu)
    else
        M[:mu] = [mu]
    end
    if cstr == :PSD
        MOM.add_constraint_nneg(M, mu)
    end
    return mu
end

function moment_variables(M, X, d::Int, cstr = :PSD)
    moment_variables(M, monomials(X,0:d), cstr)
end

