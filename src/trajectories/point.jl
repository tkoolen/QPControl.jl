struct PointTrajectory{T}
    frame::CartesianFrame3D
    traj::T
end

function (point_traj::PointTrajectory)(x, ::Val{num_derivs}) where num_derivs
    vals = point_traj.traj(x, Val(num_derivs))
    frame = point_traj.frame
    Point3D(frame, vals[1]), map(x -> FreeVector3D(frame, x), Base.tail(vals))...
end

(point_traj::PointTrajectory)(x) = first(point_traj(x, Val(0)))
