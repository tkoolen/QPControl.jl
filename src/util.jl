function val_deriv_deriv2(f, x)
    inner = ForwardDiff.Dual(x, one(x))
    outer = ForwardDiff.Dual(inner, one(x))
    ydual = f(outer)
    yddual = ForwardDiff.extract_derivative(Void, ydual)
    ForwardDiff.value(ForwardDiff.value(ydual)), ForwardDiff.value(yddual), ForwardDiff.extract_derivative(Void, yddual)
end
