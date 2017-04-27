immutable PDGains{T}
    k::T
    d::T
end

pd(gains::PDGains, e, ė) = gains.k * e + gains.d * ė
pd(gains::PDGains, x, xdes, ẋ, ẋdes) = pd(gains, xdes - x, ẋdes - ẋ)
