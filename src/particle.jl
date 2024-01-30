"Movable ball in N dimensions"
mutable struct Particle{N,T}
    x::SVector{N,T}    # position of centre
    v::SVector{N,T}    # velocity
    r::T               # radius
    m::T               # mass
    c::Bool            # collision marker
    # time_to_collision::Vector{Float64}
    # nn::Int            # next particle-particle collision
end

# Particle(x, v, r) = Particle(x, v, r, one(r), 0)
Particle(x, v, r) = Particle(x, v, r, one(r), :false)

centre(p::Particle) = p.x

"Normal vector at point x on sphere"
normal(p::Particle, x) = normalize(x - centre(p))