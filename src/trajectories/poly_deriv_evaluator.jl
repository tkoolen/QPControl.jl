# PolyDerivEvaluator: used to compute the derivatives once instead of every call
struct PolyDerivEvaluator{F, D<:Tuple{Vararg{Polynomial}}}
    f::F
    derivs::D
end

@inline function PolyDerivEvaluator(p::Polynomial, ::Val{num_derivs}) where num_derivs
    PolyDerivEvaluator(p, make_poly_derivs(p, Val(num_derivs)))
end

make_poly_derivs(p::Polynomial, ::Val{0}) = ()
@inline function make_poly_derivs(p::Polynomial, ::Val{num_derivs}) where num_derivs
    p′ = SUP.derivative(p)
    (p′, make_poly_derivs(p′, Val(num_derivs - 1))...)
end

@generated function (pde::PolyDerivEvaluator)(θ, ::Val{num_derivs}) where num_derivs
    num_derivs >= 0 || return :(throw(ArgumentError("num_derivs must be nonnegative")))
    exprs = [:(pde.f(θ))]
    for i = 1 : num_derivs
        push!(exprs, :(pde.derivs[$i](θ)))
    end
    return :(tuple($(exprs...)))
end
