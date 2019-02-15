struct SE3Trajectory{A, L}
    body::CartesianFrame3D
    base::CartesianFrame3D
    frame::CartesianFrame3D
    angular::A
    linear::L
end

function (traj::SE3Trajectory)(x, ::Val{num_derivs}) where num_derivs
    num_derivs == 2 || error()
    R, ω, ωd = traj.angular(x, Val(num_derivs))
    p, ν, νd = traj.linear(x, Val(num_derivs))
    H = Transform3D(traj.frame, traj.base, R, p)
    T = Twist(traj.body, traj.base, traj.frame, ω, ν)
    Td = SpatialAcceleration(traj.body, traj.base, traj.frame, ωd, νd)
    return H, T, Td
end
