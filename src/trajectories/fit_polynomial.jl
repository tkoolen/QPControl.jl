@inline function fit_cubic(; x0, xf, y0, yd0, yf, ydf)
    num_coeffs = Val(4)
    A = hcat(
        SVector(SUP.coefficient_gradient(x0, num_coeffs, Val(0))),
        SVector(SUP.coefficient_gradient(x0, num_coeffs, Val(1))),
        SVector(SUP.coefficient_gradient(xf, num_coeffs, Val(0))),
        SVector(SUP.coefficient_gradient(xf, num_coeffs, Val(1))),
    ) |> transpose
    b = SVector(tuple(y0, yd0, yf, ydf))
    Polynomial(Tuple(A \ b))
end

@inline function fit_quintic(; x0, xf, y0, yd0, ydd0, yf, ydf, yddf)
    num_coeffs = Val(6)
    A = hcat(
        SVector(SUP.coefficient_gradient(x0, num_coeffs, Val(0))),
        SVector(SUP.coefficient_gradient(x0, num_coeffs, Val(1))),
        SVector(SUP.coefficient_gradient(x0, num_coeffs, Val(2))),
        SVector(SUP.coefficient_gradient(xf, num_coeffs, Val(0))),
        SVector(SUP.coefficient_gradient(xf, num_coeffs, Val(1))),
        SVector(SUP.coefficient_gradient(xf, num_coeffs, Val(2))),
    ) |> transpose
    b = SVector(tuple(y0, yd0, ydd0, yf, ydf, yddf))
    Polynomial(Tuple(A \ b))
end
