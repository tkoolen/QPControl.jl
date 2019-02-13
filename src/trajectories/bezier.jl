struct BezierCurve{N, T}
    points::NTuple{N, T}
end

BezierCurve(coeffs::Tuple) = BezierCurve(promote(coeffs...))
BezierCurve(coeffs...) = BezierCurve(coeffs)

@inline StaticUnivariatePolynomials.constant(b::BezierCurve) = b.points[1]

@inline (b::BezierCurve{1})(t) = constant(b)
@inline function (b::BezierCurve)(t)
    b1 = BezierCurve(reverse(tail(reverse(b.points))))
    b2 = BezierCurve(tail(b.points))
    (oneunit(t) - t) * b1(t) + t * b2(t)
end


@inline function StaticUnivariatePolynomials.derivative(b::BezierCurve{N}) where {N}
    # https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html
    points = b.points
    n = N - 1
    BezierCurve(ntuple(i -> n * (points[i + 1] - points[i]), n))
end

"""
```math
\\int_{0}^{t} b\\left(\\tau\\right) e^{c \\tau} d\\tau
```

Found using integration by parts.
"""
@inline function exponential_integral(b::BezierCurve{N}, c, t; inv_c = inv(c), exp_c_t=exp(c * t)) where N
    # TODO: move to StaticUnivariatePolynomials?
    if N === 1
        return inv_c * constant(b) * (exp_c_t - 1)
    else
        return inv_c * (b(t) * exp_c_t - constant(b) - exponential_integral(derivative(b), c, t; inv_c=inv_c, exp_c_t=exp_c_t))
    end
end

"""
```math
\\int_{0}^{1} b\\left(\\tau\\right) e^{c \\tau} d\\tau
```
"""
@inline function exponential_integral(b::BezierCurve{N}, c; inv_c=inv(c), exp_c=exp(c)) where N
    pn = b.points[N]
    if N === 1
        return inv_c * constant(b) * (exp_c - 1)
    else
        return inv_c * (pn * exp_c - constant(b) - exponential_integral(derivative(b), c, inv_c=inv_c, exp_c=exp_c))
    end
end
