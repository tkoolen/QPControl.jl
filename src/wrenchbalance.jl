struct WrenchBalance{N}
    momentummatrix::MomentumMatrix{Matrix{Float64}}
    externalwrenches::Dict{RigidBody{Float64}, Wrench{Float64}}
    contactsettings::Vector{ContactSettings{N}}
    Ȧv_angular::Vector{Float64}
    Ȧv_linear::Vector{Float64}
end

function update!(balance::WrenchBalance, state::MechanismState)
    A = balance.momentummatrix
    centroidalframe = A.frame
    com = center_of_mass(state)
    worldframe = com.frame
    centroidal_to_world = Transform3D(centroidalframe, worldframe, com.v)
    world_to_centroidal = inv(centroidal_to_world)
    Ȧv = transform(momentum_rate_bias(state), world_to_centroidal)
    bias = Ȧv
    for wrench in values(externalwrenches)
        to_world = transform_to_root(state, wrench.frame)
        bias -= transform(wrench, world_to_centroidal * to_world)
    end
    balance.Ȧv_angular .= angular(Ȧv)
    balance.Ȧv_linear .= linear(Ȧv)
    nothing
end

function addconstraint!()
end

function constraintfunction(balance::WrenchBalance, v̇::Vector{Variable}, ρ::Vector{Variable})

end
