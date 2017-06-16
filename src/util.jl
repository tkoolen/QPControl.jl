function val_deriv_deriv2(f, x)
    ydual = f(ForwardDiff.Dual(ForwardDiff.Dual(x, one(x)), one(x)))
    yddual = ForwardDiff.extract_derivative(ydual)
    ForwardDiff.value(ForwardDiff.value(ydual)), ForwardDiff.value(yddual), ForwardDiff.extract_derivative(yddual)
end
