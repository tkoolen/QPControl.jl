struct Constant{T}
    value::T
end

make_zero_derivs(x, ::Val{0}) = ()
make_zero_derivs(x, ::Val{N}) where {N} = (deriv = zero_tangent(x); (deriv, make_zero_derivs(deriv, Val(N - 1))...))

function (traj::Constant{T})(x, ::Val{num_derivs}) where {T, num_derivs}
    traj.value, make_zero_derivs(traj.value, Val(num_derivs))...
end

# Single-argument call overload just returns the value:
(traj::Constant)(x) = first(traj(x, Val(0)))
