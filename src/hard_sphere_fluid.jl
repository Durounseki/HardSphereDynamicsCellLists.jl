"""
Hard sphere fluid in N dimensions
"""
mutable struct HardSphereFluid{N,T} <: AbstractFluid{N,T}
	box::RectangularBox{N,T}
	particles::Vector{Particle{N,T}}
end

function HardSphereFluid{N,T}(box::RectangularBox{N,T}, num_particles, r) where {N,T} #We still need a better solution so that we can use box

	particles = [Particle(zero(SVector{N,T}), zero(SVector{N,T}), r) for i in 1:num_particles]

	return HardSphereFluid{N,T}(box, particles)

end
