abstract type AbstractFluid{N,T} end

"""
Non-interacting particles fluid in N dimensions
"""
mutable struct NonInteractingFluid{N,T} <: AbstractFluid{N,T}
	box::RectangularBox{N,T}
	particles::Vector{Particle{N,T}}
end

function NonInteractingFluid{N,T}(box::RectangularBox{N,T}, num_particles, r) where {N,T}

	particles = [Particle(zero(SVector{N,T}), zero(SVector{N,T}), r) for i in 1:num_particles]

	return NonInteractingFluid{N,T}(box, particles)

end