type Maybe{T}
    valid::Bool
    data::T
    Maybe() = new(false)
    Maybe(data::T) = new(true, data)
end

Maybe{T}(data::T) = Maybe{T}(data)
isvalid(m::Maybe) = m.valid
invalidate!(m::Maybe) = (m.valid = false; nothing)
setdata!(m::Maybe, data) = (m.data = data; m.valid = true; nothing)
copydata!(m::Maybe, data) = (m.data .= data; m.valid = true; nothing)
Base.get(m::Maybe) = (isvalid(m) || error("invalid"); m.data)
