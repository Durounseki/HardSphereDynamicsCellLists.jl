function initialState(lower::SVector{N,T}, upper::SVector{N,T}, Temp::T, cell_L::T, n::Int, r::T; interaction = :none,
    conduction_alg = :isolated, radiation_alg = :isolated, τ = 0.0, ν = 0.0, Q = 0.0, kwargs...) where {N,T}

    #Create box
    box = RectangularBox(lower,upper,Temp,cell_L)
    #Define flow type, collision type and thermostat algorithms
    collision_type = ElasticCollision()
    flow_type = FreeFlow()
    conduction_thermostat = ConductionThermostat(Temp,conduction_alg)
    radiation_thermostat = RadiationThermostat(Temp, radiation_alg,τ=τ,ν=ν,Q=Q)
    #Set up simulation
    if interaction == :hard_sphere

        #Create fluid
        fluid = HardSphereFluid{N,T}(box,n,r)
        #Initialize the particle's states and the cell lists
        initial_condition!(fluid)
        event_handler = AllToAll(fluid, flow_type)

    else
        
        #Create fluid
        fluid = NonInteractingFluid{N,T}(box,n,r)
        #Initialize the particle's states and the cell lists
        initial_condition!(fluid)
        event_handler = WallToAll(fluid, flow_type)

    end

    simulation = FluidSimulation(fluid, event_handler, flow_type, collision_type, conduction_thermostat, radiation_thermostat);

    return simulation

end

function isochoric(lower::SVector{N,T}, upper::SVector{N,T}, Temp₀::T, Temp₁::T, cell_L::T, n::Int, r::T, final_time::T; δt = 0.01, steps=1, kwargs...) where {N,T}
    
    #Set initial state
    simulation = initialState(lower, upper, Temp₀, cell_L, n, r; kwargs...)
    
    @unpack fluid, conduction_thermostat, radiation_thermostat = simulation
    @unpack particles, box = fluid
    # @unpack pair_collision_times = event_handler
    
    #Save the state of the system, {particles,container}, at time t=0
    states = [deepcopy(particles)]
    times = [0.0]
    container = [deepcopy(box)]
    # collisions = [deepcopy(pair_collision_times)]

    conduction_thermostat.Temp = Temp₁
    radiation_thermostat.Temp = Temp₁

    for t in δt:steps*δt:final_time
        #Evolve system for steps*δt
        evolve!(simulation, δt, steps*δt)
        #Record last state
        push!(states,deepcopy(particles))
        push!(times,t)
        push!(container,deepcopy(box))
        # push!(collisions,deepcopy(pair_collision_times))
    end

    return simulation, states, container, times

end

function volume_expansion(lower::SVector{N,T}, upper::SVector{N,T}, Temp::T, cell_L::T, n::Int, r::T, ΔL, expansion_rate; δt = 0.01, steps=1, wall_index = 1, kwargs...) where {N,T} ;

    #Set initial state
    simulation = initialState(lower, upper, Temp, cell_L, n, r; kwargs...)

    @unpack fluid = simulation
    @unpack box, particles = fluid

    states = [deepcopy(particles)]
    times = [0.0]
    container = [deepcopy(box)]
    
    wall_displacement = expansion_rate * δt
    final_time = ΔL / expansion_rate

    moving_wall = box.walls[wall_index]

    #Add cells to the cell grid to cover the final length
    addCells(box, ΔL, wall_index = wall_index)

    for t in δt:δt:final_time
        # println(box.lower,box.upper)
        #Move the wall (for the moment try only expansion)
        moving_wall.p += wall_displacement * moving_wall.n
        # println(box.lower,box.upper)
        #update bounds
        if sum(moving_wall.n) < 0
            box.lower += wall_displacement * moving_wall.n
        else
            box.upper += wall_displacement * moving_wall.n
        end
        # println(box.lower,box.upper)

        #Calculate new state
        evolve!(simulation, δt, δt)
        push!(states, deepcopy(particles))
        push!(times, t)
        push!(container, deepcopy(box))
    end

    return simulation, states, container, times

end

function addCells(cells, L_low, L_up, L_wall, ΔL, n, cell_L)

    for i in Int(floor(L_low/cell_L)):Int(floor(L_up/cell_L))
        for j in 1:Int(ceil(ΔL/cell_L))
            cells[CartesianIndex(Int(floor(L_wall/cell_L)),i)+CartesianIndex(Int.(j*n)...)]=[]
        end
    end

end

function addCells(box,ΔL;wall_index = 1)
   
    @unpack cells, lower, upper, walls, cell_L = box
    direction = findall(x -> abs(x) > 0, walls[wall_index].n)[1] 

    return addCells(cells, lower[direction], upper[direction], walls[wall_index].p[direction], ΔL, walls[wall_index].n, cell_L)

end


function adiabatic_expansion(lower::SVector{N,T}, upper::SVector{N,T}, Temp₀::T, Temp₁, cell_L::T, n::Int, r::T, ΔL, expansion_rate;
    δt = 0.01, steps=100, wall_index = 1, radiation_alg = :v_scaling, kwargs...) where {N,T} ;

    #Set initial state
    simulation = initialState(lower, upper, Temp₀, cell_L, n, r; radiation_alg = radiation_alg, kwargs...)

    @unpack fluid, conduction_thermostat, radiation_thermostat = simulation
    @unpack box, particles = fluid

    states = [deepcopy(particles)]
    times = [0.0]
    container = [deepcopy(box)]
    
    wall_displacement = expansion_rate * δt
    final_time = (ΔL / expansion_rate)*steps

    moving_wall = box.walls[wall_index]

    #Add cells to the cell grid to cover the final length
    addCells(box, ΔL, wall_index = wall_index)

    #Adiabatic condition
    C₀ = adiabaticConstant(box.lower,box.upper,Temp₀)
    Temp₁ = Temp₀

    for t in δt:steps*δt:final_time

        #Move the wall (for the moment try only expansion)
        moving_wall.p += wall_displacement * moving_wall.n

        #update bounds
        if sum(moving_wall.n) < 0
            box.lower += wall_displacement * moving_wall.n
        else
            box.upper += wall_displacement * moving_wall.n
        end

        # Let's use velocity rescaling for the moment
        conduction_thermostat.Temp = Temp₁
        radiation_thermostat.Temp = Temp₁

        #For the moment we consider only the v scaling thermostat
        #Calculate state at new volume and equilibrate for steps number of timesteps
        for i in 1:steps
            evolve!(simulation, δt, δt)
            push!(states, deepcopy(particles))
            push!(times, t+(i-1)*δt)
            push!(container, deepcopy(box))
        end
        #Calculate new temperature
        Temp₁ = adiabaticCondition(box.lower,box.upper,C₀)
    end

    return simulation, states, container, times

end

function adiabaticCondition(lower::SVector{N,T},upper::SVector{N,T},C₀::Float64) where {N,T}
    Vol = prod(upper-lower)
    γ = (1.0*N+2)/N
    return C₀/Vol^(γ-1)
end

function adiabaticConstant(lower::SVector{N,T},upper::SVector{N,T},Temp::Float64) where {N,T}
    γ = (1.0*N+2)/N
    Vol = prod(upper-lower)
    return Temp*Vol^(γ-1)
end