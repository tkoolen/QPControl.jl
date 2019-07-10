struct PointTrajectory{T}
    frame::CartesianFrame3D
    trajectory::T
end

function (traj::PointTrajectory)(x, ::Val{num_derivs}) where num_derivs
    vals = traj.trajectory(x, Val(num_derivs))
    frame = traj.frame
    Point3D(frame, vals[1]), map(x -> FreeVector3D(frame, x), Base.tail(vals))...
end

(traj::PointTrajectory)(x) = Point3D(traj.frame, traj.trajectory(x))


struct FreeVectorTrajectory{T}
    frame::CartesianFrame3D
    trajectory::T
end

function (traj::FreeVectorTrajectory)(x, ::Val{num_derivs}) where num_derivs
    vals = traj.trajectory(x, Val(num_derivs))
    frame = traj.frame
    map(x -> FreeVector3D(frame, x), vals)
end

(traj::FreeVectorTrajectory)(x) = FreeVector3D(traj.frame, traj.trajectory(x))
