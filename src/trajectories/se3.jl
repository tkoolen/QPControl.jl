struct SE3Trajectory{A, L}
    body::CartesianFrame3D
    base::CartesianFrame3D
    angular::A # rotation from body to base, derivatives expressed in base
    linear::L # translation from body to base, derivatives expressed in bsae
end

function (traj::SE3Trajectory)(x, ::Val{num_derivs}) where num_derivs
    num_derivs == 2 || error()
    R, ω_base, ωd_base = traj.angular(x, Val(num_derivs))
    p, pd, pdd = traj.linear(x, Val(num_derivs))
    Rinv = inv(R)

    H = Transform3D(traj.body, traj.base, R, p)

    ω_body = Rinv * ω_base
    ν_body = Rinv * pd
    T = Twist(traj.body, traj.base, traj.body, ω_body, ν_body)

    ωd_body = Rinv * ωd_base
    νd_body = Rinv * pdd + ω_body × ν_body

    Td = SpatialAcceleration(traj.body, traj.base, traj.body, ωd_body, νd_body)
    return H, T, Td
end
