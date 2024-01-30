abstract type AbstractThermostat end

"""

Assumes the gas container, the box, is in contact with a thermal bath and each particle thermalizes every time it collides with a wall.

`Temp`: The temperature of the heat bath.
`algorithm`: Determines the resulting speed of the particle after a collision.
                    `:isolated` No change of the particle's speed after a collision.
                    `:mb_dist_sampling` Samples a new speed from the Maxwell-Boltzmann distribution at temperature Temp.
                    `:constant_sampling` Enforces the average kinetic energy to remain constant.
                    `:mb_normal` samples only the normal component of the velocity to the wall. 

"""
mutable struct ConductionThermostat{Float64,Symbol} <: AbstractThermostat
    Temp::Float64
    algorithm::Symbol
end

"""

Dictionary of symbols representing the different conduction thermostat algorithms.

"""
ConductionAlgorithms = Dict{Int,Symbol}([(1,:mb_dist_sampling),(2,:constant_sampling),(3,:mb_normal)])

"""

Builds the collision thermostat using the dictionary of symbols.
`Temp`: The temperature of the heat bath
`i`:    1 -> :mb_dist_sampling
        2 -> :constant_sampling
        3 -> :mb_normal

"""

ConductionThermostat{Float64,Int}(Temp,i) = ConductionThermostat(Temp,ConductionAlgorithms[i])

"""

Builds the conduction thermostat using the symbol for the algorithm directly.
`Temp`: The temperature of the heat bath
`algorithm`: Determines the resulting speed of the particle after a collision.
                    `:isolated` No change of the particle's speed after a collision.
                    `:mb_dist_sampling` Samples a new speed from the Maxwell-Boltzmann distribution at temperature Temp.
                    `:constant_sampling` Enforces the average kinetic energy to remain constant.
                    `:mb_normal` samples only the normal component of the velocity to the wall. 

"""

ConductionThermostat{Float64,Symbol}(Temp,algorithm) = ConductionThermostat(Temp,algorithm)

"""

Builds the collision thermostat with the default algorithm for an isolated system.

"""

ConductionThermostat() = ConductionThermostat(0.0,:isolated)

"""

Assumes the gas container, the box, is in contact with a thermal bath and each particle thermalizes through radiation.

`Temp`: The temperature of the heat bath.
`algorithm`: Determines the new velocities of the ensemble of particles
                    `:isolated` No change of the particle's speed.
                    `:v_scaling` The speeds of all the particles are rescaled by a factor λ₀ = sqrt( K(T_bath) / K(T) ), where K(T) = (N / 2) T, is the kinetic energy.
                    `:berendsen` Berendsen thermostat. The velocity is rescaled over a period of time τ, λ₁ = sqrt( 1 + (δt/τ)*( λ₀ - 1 ) )
                    `:stochastic_v_scaling` Similar to `:vscaling` but with K(T_bath) = K_b drawn from the gamma distribution.
                    `:b_d_p` Bussi-Donadio-Parrinello thermostat, similar to the Berendsen method, but with K_b drawn from a gamma distribution.
                    `:andersen` Andersen thermostat. A new speed is sampled from the Maxwell-Boltzmann distribution through a Poisson process with rate ν.
                    `:nose_hoover` Nosè-Hoover thermostat. The hamiltonian of the system includes the bath explicitly in the hamiltonian. The bath has virtual mass Q.

"""

mutable struct RadiationThermostat{Float64,Symbol} <: AbstractThermostat
    Temp::Float64
    algorithm::Symbol
    τ::Float64
    ν::Float64
    Q::Float64
    function RadiationThermostat(Temp,algorithm; kwargs...)
        #Set default values for the algorithm parameters
        this = new{Float64,Symbol}(Temp,algorithm, 0.0, 0.0, 0.0)
        for (key, value) in Iterators.pairs(kwargs)
            setfield!(this,key,value)
        end
        # if haskey(kwargs,:τ )
        #     this.τ = kwargs[:τ]
        # else
        #     this.τ = 0.0
        # end
        # if haskey(kwargs,:ν)
        #     this.ν = kwargs[:ν]
        # else
        #     this.ν = 0.0
        # end
        # if haskey(kwargs,:Q)
        #     this.Q = kwargs[:Q]
        # else
        #     this.Q = 0.0
        # end
        return this
    end
end

"""

Dictionary of symbols representing the different radiation thermostat algorithms.

"""

RadiationAlgorithms = Dict{Int,Symbol}([(1,:v_scaling), (2,:berendsen), (3,:stochastic_v_scaling), (4,:b_d_p), (5,:andersen), (6,:nose_hoover)])

"""

Builds the collision thermostat using the dictionary of symbols.
`Temp`: The temperature of the heat bath
`i`:    1 -> :v_scaling
        2 -> :berendsen, parameters τ -> kw = τ 
        3 -> :stochastic_v_scaling
        4 -> :b_d_p
        5 -> :andersen
        6 -> :nose_hoover

"""

# RadiationThermostat{Float64,Int}(Temp,i; kwargs...) = RadiationThermostat(;Temp=Temp, algorithm = RadiationAlgorithms[i], kwargs...)
RadiationThermostat{Float64,Int}(Temp,i; kwargs...) = RadiationThermostat(Temp, RadiationAlgorithms[i], kwargs...)

# """

# Builds the radiation thermostat using the symbol for the algorithm directly.
# `Temp`: The temperature of the heat bath
# `algorithm`: Determines the new velocities of the ensemble of particles
#                     `:isolated` No change of the particle's speed.
#                     `:v_scaling` The speeds of all the particles are rescaled by a factor λ₀ = sqrt( K(T_bath) / K(T) ), where K(T) = (N / 2) T, is the kinetic energy.
#                     `:berendsen` Berendsen thermostat. The velocity is rescaled over a period of time τ, λ₁ = sqrt( 1 + (δt/τ)*( λ₀ - 1 ) )
#                     `:stochastic_v_scaling` Similar to `:vscaling` but with K(T_bath) = K_b drawn from the gamma distribution.
#                     `:b_d_p` Bussi-Donadio-Parrinello thermostat, similar to the Berendsen method, but with K_b drawn from a gamma distribution.
#                     `:andersen` Andersen thermostat. A new speed is sampled from the Maxwell-Boltzmann distribution through a Poisson process with rate ν.
#                     `:nose_hoover` Nosè-Hoover thermostat. The hamiltonian of the system includes the bath explicitly in the hamiltonian. The bath has virtual mass Q.

# """

# # RadiationThermostat{Float64,Symbol}(Temp,algorithm; kwargs...) = RadiationThermostat(;Temp=Temp ,algorithm=algorithm, kwargs...)
# RadiationThermostat{Float64,Symbol}(Temp,algorithm; kwargs...) = RadiationThermostat(Temp=Temp, algorithm=algorithm, kwargs...)

# # """

# # Builds the radiation thermostat with the default algorithm for an isolated system.

# # """

# # RadiationThermostat() = RadiationThermostat(0.0,:isolated,0.0,0.0,0.0)

"""

Thermalizes the particles that collided with a wall.
`particles`: array of the particles in the system.
`thermostat`: CollisionThermostat at temperature `Temp`

"""
function thermalize!(particles::Vector{Particle{N,T}}, thermostat::ConductionThermostat) where {N,T}
    
    @unpack algorithm, Temp = thermostat

    #Find the particles that collided
    reflected_particles = particles[findall( p -> p.c, particles)]

    #Sample the new speeds for the particles that collided
    if algorithm == :mb_dist_sampling
        #Sample new speed from the maxwell boltzmann distribution
        new_v = [norm(SA[[MB_dist(Temp) for i in 1:N]...]) for p in reflected_particles]

    elseif algorithm == :constant_sampling
        #Average kinetic energy = (N/2) * Temperature, with kB = 1
        new_v = [sqrt( ( N / p.m ) * Temp ) for p in reflected_particles]
    end
    
    for i in eachindex(reflected_particles)
        reflected_particles[i].v *= new_v[i] / norm(reflected_particles[i].v)
        reflected_particles[i].c = :false
    end

end

"""

Thermalizes all the particles through radiation.
`particles`: array of the particles in the system.
`thermostat`: RadiationThermostat at temperature `Temp`
`δt` time step

"""

function thermalize!(particles::Vector{Particle{N,T}}, thermostat::RadiationThermostat, δt) where {N,T}

    @unpack Temp, algorithm, τ, ν, Q = thermostat
    
	T_state = sum([ p.m * normsq(p.v) for p in particles]) / (2*length(particles))

    if algorithm == :v_scaling

        for p in particles
            p.v *= sqrt( Temp / T_state )
        end

	elseif algorithm == :berendsen

        for p in particles
            p.v *= sqrt( 1 + ( δt / τ ) * ( Temp / T_state - 1 )  )
        end

    elseif algorithm == :stochastic_v_scaling

        T_bath = sum([p.m * normsq(SA[MB_dist(Temp), MB_dist(Temp)]) for p in particles]) / (2 * length(particles))

        # T_bath = ( 2.0 / N ) * rand(Gamma( N / 2.0 , Temp ) )#Sample from gamma with β = (1/radiation_thermostat.Temp)

        for p in particles
            p.v *= sqrt( T_bath / T_state )
        end

    elseif algorithm == :b_d_p

        #Calculate the Kinetic energy
        K_state = ( N / 2.0 ) * T_state
        K_bath = ( N / 2.0 ) * Temp
        Rs = [randn() for i in 1:N]

        α = sqrt( exp(-δt / τ) + ( K_bath / (N * K_state) ) * (1-exp(-δt / τ)) * sum(Rs .^ 2) + 2 * exp(-δt / (2*τ)) * sqrt( ( K_bath / (N * K_state) ) * (1-exp(-δt / τ)) ) * Rs[1] )
        # if (Rs[1] + sqrt(2 * K_bath / K_state * exp(-δt / τ) / (1-exp(-δt / τ))) ) < 0
        #     α *= -1
        # end
        
        #Conserved energy quantity
        # H = K * ( 1 - α ^2 )

        for p in particles
            p.v *= α
        end

    elseif algorithm == :andersen

        for p in particles
            
            if rand() < ν*δt
                p.v = SA[[MB_dist(Temp) for i in 1:N]...]
                # p.v *=  norm(SA[[MB_dist(Temp) for i in 1:N]...]) / norm(p.v)
                # p.c = :true
            end
        end

    end

end

