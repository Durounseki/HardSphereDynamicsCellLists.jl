abstract type AbstractCollisionDynamics end

struct ElasticCollision <: AbstractCollisionDynamics end


"""

Elastic collision of ball with Wall. The ball is assumed to be touching the Wall. The speed of the ball after the collision is modified according to a thermalization
algorithm.

`b`: particle in N dimensions with position x and velocity v.
`Π`: fixed wall.
`thermostat`: defines the type of thermalization algorithm after a particle collides with a wall.


"""

function collide!(b::Particle{N,T}, Π::Wall{N,T}, ::ElasticCollision) where {N,T}
    
    v = b.v
    n = Π.n

    b.v = v - 2 * (v⋅n) * n

    b.c = :true

    b.col_time = Inf
    b.col_pair = -1
    b.col_type = :none

end


"Assumes b1 and b2 are touching"
function collide!(b1::Particle, b2::Particle, ::ElasticCollision)
    Δx = b1.x - b2.x

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

    b1.col_time = Inf
    b1.col_pair = -1
    b1.col_type = :none
    b2.col_time = Inf
    b2.col_pair = -1
    b2.col_type = :none

end

#The maxwell boltzmann distribution in one dimension
function MB_dist(T;m=1,kB=1)
    return sqrt(kB*T/m)randn()
end