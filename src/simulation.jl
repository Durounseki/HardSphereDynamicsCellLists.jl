mutable struct FluidSimulation{
    N, T,
    F<:AbstractFluid{N,T},
    H<:AbstractEventHandler,
    L<:AbstractFlowDynamics,
    C<:AbstractCollisionDynamics,
	O<:AbstractThermostat,
	R<:AbstractThermostat
	}

    fluid::F
    event_handler::H
    flow_dynamics::L
    collision_dynamics::C
	conduction_thermostat::O
	radiation_thermostat::R
	current_time::T

end

function FluidSimulation(fluid::AbstractFluid{N,T}, event_handler, flow_dynamics, collision_dynamics, conduction_thermostat, radiation_thermostat) where {N,T}
    FluidSimulation(fluid, event_handler, flow_dynamics, collision_dynamics, conduction_thermostat, radiation_thermostat, zero(T))
end

# mutable struct NonInteractingSimulation{
# 	N,T,
# 	F<:NonInteractingFluid{N,T},
#     H<:AbstractEventHandler,
#     L<:AbstractFlowDynamics,
#     C<:AbstractCollisionDynamics,
# 	O<:AbstractThermostat,
# 	R<:AbstractThermostat
# 	}
# 	fluid::F
#     event_handler::H
#     flow_dynamics::L
#     collision_dynamics::C
# 	conduction_thermostat::O
# 	radiation_thermostat::R
# 	current_time::T
# end

# function NonInteractingSimulation(fluid::NonInteractingFluid{N,T}, event_handler, flow_dynamics, collision_dynamics, conduction_thermostat, radiation_thermostat) where {N,T}
#     NonInteractingSimulation(fluid, event_handler, flow_dynamics, collision_dynamics, conduction_thermostat, radiation_thermostat, zero(T))
# end

function normsq(v)
    return sum(abs2, v)
end

function overlap(b1::Particle, b2::Particle)
    normsq(b1.x - b2.x) < (b1.r + b2.r)^2
end

function overlap(b::Particle, particles::Vector{Particle{N,T}}) where {N,T}
    any(overlap.(Ref(b), particles))
end


"""
Generate initial condition for particles.

Uses random sequential deposition:
place a disc at a time and check that it doesn't overlap with any previously placed disc.

Generates allowed ball positions and uniform velocities.
"""
function initial_condition!(particles::Vector{Particle{N,T}}, table, fluid_type; equilibrium = :true,
		lower=table.lower, upper=table.upper) where {N,T}

	## TODO: Check that initial condition is OK with respect to planes

    U = Uniform.(lower .+ particles[1].r, upper .- particles[1].r)
    particles[1].x = rand.(U)

    for i = 2:length(particles)
        U = Uniform.(lower .+ particles[i].r, upper .- particles[i].r)
        particles[i].x = rand.(U)

		if fluid_type <: HardSphereFluid
			count = 0

			while overlap(particles[i], particles[1:i-1])
				particles[i].x = rand.(U)
				count += 1

				count > 10^5 && error("Unable to place disc $i")
			end
		end

    end

	if equilibrium
		
		#Sample the components of the velocities from the Maxwell-Boltzmann distribution
    	for i = 1:length(particles)
        	particles[i].v = SA[MB_dist(table.Temp),MB_dist(table.Temp)]
    	end

	else
		
		for i = 1:length(particles)
        	particles[i].v = sqrt(N*table.Temp/particles[i].m)*normalize(@SVector randn(N))
		end
    
	end

    #Generate initial cell lists
    CellList([p.x for p in particles],table)

end

#Implementation of initial equilibrium condition missing
initial_condition!(fluid; kw...) = initial_condition!(fluid.particles, fluid.box, typeof(fluid); kw...)

"Carry out collision assuming already at moment of collision"
function collide!(fluid::AbstractFluid, event_handler, collision_dynamics)

	@unpack particles, box = fluid
	@unpack partner1, partner2, collision_type = event_handler

	if collision_type == :wall_collision
		collide!(particles[partner1], box.walls[partner2], collision_dynamics)

	elseif collision_type == :disc_collision
		collide!(particles[partner1], particles[partner2], collision_dynamics)

	else
		error("No collision")
	end

end

# #Non Interacting gas collisions
# "Carry out collision assuming already at moment of collision"
# function collide!(fluid::NonInteractingFluid, event_handler, collision_dynamics)

# 	@unpack particles, box = fluid
# 	@unpack partner1, partner2, collision_type = event_handler

# 	collide!(particles[partner1], box.walls[partner2], collision_dynamics)
# end



"""
	evolve!(fluid::HardSphereFluid, num_collisions::Integer)

Discrete-time evolution: Evolve fluid for `num_collisions` collisions.
Both sphere--sphere and sphere--wall collisions are counted as collisions.

Returns post-collision states, times at which collisions occur, and collision types.
"""
function evolve!(simulation::FluidSimulation{N,T}, num_collisions::Integer) where {N,T}

	@unpack fluid, flow_dynamics, collision_dynamics, event_handler = simulation

	states = [deepcopy(fluid.particles)]
	times = [0.0]
	collision_types = [:none]


    for i in 1:num_collisions

		flow!(fluid.particles, fluid.box, (event_handler.next_collision_time - simulation.current_time), flow_dynamics)
		simulation.current_time = event_handler.next_collision_time

		push!(collision_types, event_handler.collision_type)

		collide!(fluid, event_handler, collision_dynamics)
		find_collision!(event_handler, fluid, flow_dynamics)

		push!(states, deepcopy(fluid.particles))
        push!(times, simulation.current_time)

    end

    return states, times, collision_types

end


"""
	flow!(fluid::HardSphereFluid, t)

Update hard sphere fluid by flowing for a time `t`, taking
account of possible collisions during that time.
"""

function flow!(simulation::FluidSimulation, t)

	@unpack fluid, flow_dynamics, collision_dynamics, event_handler, conduction_thermostat = simulation
	# @unpack particles, table = fluid

	time_to_next_collision = event_handler.next_collision_time - simulation.current_time

	while t > time_to_next_collision
		t -= time_to_next_collision

		flow!(fluid.particles, fluid.box, time_to_next_collision, flow_dynamics)
		simulation.current_time += time_to_next_collision

		collide!(fluid, event_handler, collision_dynamics)
		find_collision!(event_handler, fluid, flow_dynamics)

		time_to_next_collision = event_handler.next_collision_time - simulation.current_time

	end

	flow!(fluid.particles, fluid.box, t, flow_dynamics)
	simulation.current_time += t

end


"""
	evolve!(fluid::HardSphereFluid, times)

Time evolution, calculating positions and velocities at given times.
"""
function evolve!(simulation::FluidSimulation{N,T}, times) where {N,T}

	@unpack fluid, flow_dynamics, collision_dynamics, event_handler, conduction_thermostat, radiation_thermostat = simulation

	# states = [deepcopy(fluid.particles)]
	# ts = [0.0]
	# cell_lists = [deepcopy(fluid.box.cells)] 

	current_t = 0.0

	for t in times
		flow!(simulation, t - current_t)

		#Thermalize the system through conduction
		if conduction_thermostat.algorithm != :isolated
			thermalize!(fluid.particles, conduction_thermostat)
		end
		#Thermalize the system through radiation
		if radiation_thermostat.algorithm != :isolated
			thermalize!(fluid.particles, radiation_thermostat, Float64(times.step))
		end

        # push!(states, deepcopy(fluid.particles))
		# push!(ts, t)
		# push!(cell_lists, deepcopy(fluid.box.cells))

		current_t = t
    end

    # return states, ts, cell_lists
end


"""
	evolve!(fluid::HardSphereFluid, δt, final_time)

Calculate positions and velocities at times up to `final_time`,
spaced by `δt`.

"""
evolve!(simulation::FluidSimulation, δt, final_time) = evolve!(simulation, δt:δt:final_time)

