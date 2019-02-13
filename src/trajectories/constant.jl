struct Constant{T}
    value::T
end

function (traj::Constant{T})(x, ::Val{num_derivs}) where {T, num_derivs}
    traj.value, ntuple(_ -> zero(tangent_type(T)), Val(num_derivs))...
end

# Single-argument call overload just returns the value:
(traj::Constant)(x) = first(traj(x, Val(0)))
