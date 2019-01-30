struct Constant{T}
    value::T
end

function (traj::Constant{T})(x, ::Val{num_derivs}=Val(0)) where {T, num_derivs}
    traj.value, ntuple(_ -> zero(tangent_type(T)), Val(num_derivs))...
end
