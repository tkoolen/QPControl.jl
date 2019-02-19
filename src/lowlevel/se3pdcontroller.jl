struct SE3PDController{T, W, G<:SE3PDGains}
    base::BodyID
    body::BodyID
    trajectory::Base.RefValue{T}
    weight::W
    gains::Base.RefValue{G}
end

function SE3PDController(base::BodyID, body::BodyID, trajectory, weight, gains::SE3PDGains)
    SE3PDController(base, body, Ref(trajectory), weight, Ref(gains))
end

function (controller::SE3PDController)(t, state::MechanismState)
    Href, Tref, Tdref = controller.trajectory[](t, Val(2))
    H = relative_transform(state, Href.from, Href.to)
    T = transform(relative_twist(state, controller.body, controller.base), inv(transform_to_root(state, Tref.frame)))
    return Tdref + pd(controller.gains[], H, Href, T, Tref)
end
