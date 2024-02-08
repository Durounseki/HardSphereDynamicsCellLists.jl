abstract type AbstractEventHandler end

mutable struct AllToAll{T} <: AbstractEventHandler
    next_collision_time::T
    partner1::Int
    partner2::Int
    collision_type::Symbol
	# pair_collision_times::Array{T,2}
end

AllToAll{T}() where {T} = AllToAll{T}(0, -1, -1, :none)
# AllToAll{T}(n) where {T} = AllToAll{T}(0, -1, -1, :none, [Inf for i in 1:n, j in 1:n])

function AllToAll(fluid::HardSphereFluid{N,T}, flow_type) where {N,T}
	event_handler = AllToAll{T}()
	find_collision!(event_handler, fluid, flow_type)

	return event_handler
end

mutable struct WallToAll{T} <: AbstractEventHandler
    next_collision_time::T
    partner1::Int
    partner2::Int
	collision_type::Symbol
end

WallToAll{T}() where{T} = WallToAll{T}(0,-1,-1,:wall_collision)

function WallToAll(fluid::NonInteractingFluid{N,T}, flow_type) where {N,T}
	event_handler = WallToAll{T}()
	find_collision!(event_handler, fluid, flow_type)

	return event_handler
end

# function initial_collision_times()



# end
"""
find the next collision event, wall-particle or particle-particle, for a hard-sphere fluid
"""
function find_collision(event_handler::AllToAll, particles, box, flow::AbstractFlowDynamics)

	partner1 = -1
	partner2 = -1
	min_collision_time = Inf
	collision_type = :none
	n=length(particles)

	for i in eachindex(particles), j in eachindex(box.walls)

		t = collision_time(particles[i], box.walls[j], flow)

		if t < particles[i].col_time
			particles[i].col_time = t
			particles[i].col_pair = n + j
		end

	end

	#Calculate the next colliding pair of particles
    nextCollidingPair(box.cells, particles, flow)

	for i in eachindex(particles)
		if particles[i].col_time < min_collision_time
			min_collision_time = particles[i].col_time
			partner1 = i
			if particles[i].col_pair <= n
				partner2 = particles[i].col_pair
				collision_type = :disc_collision
			else
				partner2 = particles[i].col_pair - n
				collision_type = :wall_collision
			end
		end
	end

	return partner1, partner2, collision_type, min_collision_time

end


find_collision(event_handler::AllToAll, fluid::HardSphereFluid, flow_type) = find_collision(event_handler, fluid.particles, fluid.box, flow_type)


function find_collision!(event_handler::AllToAll, fluid::HardSphereFluid, flow_type)
	partner1, partner2, collision_type, min_collision_time = find_collision(event_handler, fluid, flow_type)

	event_handler.partner1 = partner1
	event_handler.partner2 = partner2
	event_handler.collision_type = collision_type
	event_handler.next_collision_time += min_collision_time
end

"""
Find the next collision, wall-particle, for a non-interacting particle fluid
"""
function find_collision(event_handler::WallToAll, particles, box, flow::AbstractFlowDynamics)

	partner1 = -1
	partner2 = -1
	min_collision_time = Inf

	for i in eachindex(particles), j in eachindex(box.walls)

		t = collision_time(particles[i], box.walls[j], flow)

		if t < particles[i].col_time
			particles[i].col_time = t
			particles[i].col_pair = j
		end

	end

	for i in eachindex(particles)
		if particles[i].col_time < min_collision_time
			min_collision_time = particles[i].col_time
			partner1 = i
			partner2 = particles[i].col_pair
		end
	end

	return partner1, partner2, min_collision_time

end

find_collision(event_handler::WallToAll, fluid::NonInteractingFluid, flow_type) = find_collision(event_handler, fluid.particles, fluid.box, flow_type)

function find_collision!(event_handler::WallToAll, fluid::NonInteractingFluid, flow_type)
	
	partner1, partner2, min_collision_time = find_collision(event_handler, fluid, flow_type)
	event_handler.partner1 = partner1
	event_handler.partner2 = partner2
	event_handler.next_collision_time += min_collision_time

end