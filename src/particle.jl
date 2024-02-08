"Movable ball in N dimensions"
mutable struct Particle{N,T}
    x::SVector{N,T}    # position of centre
    v::SVector{N,T}    # velocity
    r::T               # radius
    m::T               # mass
    c::Bool            # collision marker
    col_time::T        # time to next collision
    col_pair::Int      # next colliding neighbor
    col_type::Symbol   # collision type
end

# Particle(x, v, r) = Particle(x, v, r, one(r), 0)
# Particle(x, v, r) = Particle(x, v, r, one(r), :false)
Particle(x, v, r) = Particle(x, v, r, one(r), :false, Inf, -1, :nothing)

centre(p::Particle) = p.x

"Normal vector at point x on sphere"
normal(p::Particle, x) = normalize(x - centre(p))