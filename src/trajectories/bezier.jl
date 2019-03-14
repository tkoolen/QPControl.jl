
struct BezierCurve{N, T} <: SUP.AbstractPolynomial{N, T}
    points::NTuple{N, T}
end

BezierCurve(coeffs::Tuple) = BezierCurve(promote(coeffs...))
BezierCurve(coeffs...) = BezierCurve(coeffs)

# Conversion to monomial basis
Polynomial(b::BezierCurve{1}) = Polynomial(b.points)
function Polynomial(b::BezierCurve{N, T}) where {N, T}
    # b(t) = (1 - t) * b1(t) + t * b2(t)
    #      = b1(t) + t * (b2(t) - b1(t))
    b1 = BezierCurve(reverse(tail(reverse(b.points))))
    b2 = BezierCurve(tail(b.points))
    Polynomial(b1) + Polynomial(0, 1) * Polynomial(b2 - b1)
end

# Utility
@inline SUP.constant(b::BezierCurve) = b.points[1]
Base.zero(::Type{BezierCurve{N, T}}) where {N, T} = BezierCurve(ntuple(_ -> zero(T), Val(N)))
Base.zero(b::BezierCurve) = zero(typeof(b))

@inline (b::BezierCurve{1})(t) = constant(b)
@inline function (b::BezierCurve)(t)
    b1 = BezierCurve(reverse(tail(reverse(b.points))))
    b2 = BezierCurve(tail(b.points))
    (oneunit(t) - t) * b1(t) + t * b2(t)
end

for op in [:+, :-]
    @eval begin
        # Two BezierCurves
        Base.$op(b1::BezierCurve{N}, b2::BezierCurve{N}) where {N} = BezierCurve(_map($op, b1.points, b2.points))

        # BezierCurves and constant
        Base.$op(b::BezierCurve, c) = BezierCurve(_map(p -> $op(p, c), b.points))
        Base.$op(c, b::BezierCurve) = BezierCurve(_map(p -> $op(c, p), b.points))
    end
end

for op in [:*, :/]
    @eval Base.$op(b::BezierCurve, c) = BezierCurve(_map(x -> $op(x, c), b.points))
end
Base.:*(c, b::BezierCurve) = BezierCurve(_map(x -> c * x, b.points))

@inline function SUP.derivative(b::BezierCurve{N}) where {N}
    # https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html
    points = b.points
    n = N - 1
    BezierCurve(ntuple(i -> n * (points[i + 1] - points[i]), n))
end


"""
```math
\\int_{0}^{1} b\\left(\\tau\\right) e^{c \\tau} d\\tau
```
"""
@inline function SUP.exponential_integral(b::BezierCurve{N}, c; inv_c=inv(c), exp_c=exp(c)) where N
    pn = b.points[N]
    if N === 1
        return inv_c * constant(b) * (exp_c - 1)
    else
        return inv_c * (pn * exp_c - constant(b) - exponential_integral(derivative(b), c, inv_c=inv_c, exp_c=exp_c))
    end
end
