abstract type AbstractEventHandler end

mutable struct AllToAll{T} <: AbstractEventHandler
    next_collision_time::T
    partner1::Int
    partner2::Int
    collision_type::Symbol
	pair_collision_times::Array{T,2}
end

AllToAll{T}(n) where {T} = AllToAll{T}(0, -1, -1, :none, [Inf for i in 1:n, j in 1:n])

function AllToAll(fluid::HardSphereFluid{N,T}, flow_type) where {N,T}
	event_handler = AllToAll{T}(length(fluid.particles))
	find_collision!(event_handler, fluid, flow_type)

	return event_handler
end

function find_collision(event_handler::AllToAll, particles, box, flow::AbstractFlowDynamics)

	partner1 = -1
	partner2 = -1
	min_collision_time = Inf
	collision_type = :none

	for i in 1:length(particles), j in 1:length(box.walls)

		t = collision_time(particles[i], box.walls[j], flow)

		if t < min_collision_time
			partner1 = i
			partner2 = j
			collision_type = :wall_collision
			min_collision_time = t
		end

	end

    #Calculate the next colliding pair of particles
    t, i, j = nextCollidingPair(box.cells, particles, flow, event_handler.pair_collision_times)

    if t < min_collision_time
        partner1 = i
        partner2 = j
        collision_type = :disc_collision
        min_collision_time = t
    end

	return partner1, partner2, collision_type, min_collision_time

end


find_collision(event_handler, fluid::HardSphereFluid, flow_type) = find_collision(event_handler, fluid.particles, fluid.box, flow_type)


function find_collision!(event_handler::AllToAll, fluid::HardSphereFluid, flow_type)
	partner1, partner2, collision_type, min_collision_time = find_collision(event_handler, fluid, flow_type)

	event_handler.partner1 = partner1
	event_handler.partner2 = partner2
	event_handler.collision_type = collision_type
	event_handler.next_collision_time += min_collision_time
end