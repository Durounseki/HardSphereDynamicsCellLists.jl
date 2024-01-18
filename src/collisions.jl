abstract type AbstractCollisionDynamics end

struct ElasticCollision <: AbstractCollisionDynamics end



"Elastic collision of ball with Wall.
The ball is assumed to be touching the Wall."

function collide!(b::Particle{N,T}, Π::Wall{N,T}, ::ElasticCollision) where {N,T}
    v = b.v
    #sample new speed from the maxwell boltzmann distribution
    new_v = sqrt(MB_dist(Π.Temp)^2+MB_dist(Π.Temp)^2)
    n = Π.n

    b.v = new_v*normalize(v - 2 * (v⋅n) * n)
end


"Assumes b1 and b2 are touching"
function collide!(b1::Particle, b2::Particle, ::ElasticCollision)
    Δx = b1.x - b2.x

    # @show norm(Δx)

    v1 = b1.v
    v2 = b2.v

    m1 = b1.m
    m2 = b2.m

    # reference: https://en.wikipedia.org/wiki/Elastic_collision#Two-dimensional_collision_with_two_moving_objects

    # Component orthogonal to line joining centers
    # is unchanged by elastic collision (since impulse is applied in direction joining sphere centers)

    n = normalize(Δx)   # vector joining sphere centers

    v1n = (v1 ⋅ n) * n  # velocity components in that direction
    v2n = (v2 ⋅ n) * n

    # treat terms m1 / (m1 + m2)  and m2 / (m1 + m2)
    # by dividing top and bottom by m1
    # This allows one of them to be Inf

    mass_ratio = (1 / ((m1 / m2) + 1))

    v1′ = v1 + 2 * mass_ratio * (v2n - v1n)  # update only components in n direction
    v2′ = v2 + 2 * (1 - mass_ratio) * (v1n - v2n)

    b1.v = v1′
    b2.v = v2′

end

#The maxwell boltzmann distribution in one dimension
function MB_dist(T;m=1,kB=1)
    return sqrt(kB*T/m)randn()
end