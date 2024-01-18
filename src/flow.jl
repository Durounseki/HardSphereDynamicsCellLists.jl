abstract type AbstractFlowDynamics end


function flow!(particles::Vector{<:Particle}, cells, cell_L, t, flow_type::AbstractFlowDynamics)
	
    for particle in particles
		flow!(particle, t, flow_type)
	end

    #calculate new cell lists
    emptyCells(cells)
    CellList([p.x for p in particles],cell_L,cells)

end


struct FreeFlow <: AbstractFlowDynamics end


function collision_time(b::Particle, Π::Wall, ::FreeFlow)
    @unpack r, x, v = b
    @unpack n, p, Temp = Π

    if v ⋅ n < 0  # travelling in wrong direction
        return Inf
    end

    return (-r - (x - p)⋅n) / (v ⋅ n)
end



function collision_time(b1::Particle, b2::Particle, ::FreeFlow)
    Δx = b1.x - b2.x
    Δv = b1.v - b2.v

    b = Δx⋅Δv

    b > 0 && return Inf  # moving away from one another


    a = normsq(Δv)
    c = normsq(Δx) - (b1.r + b2.r)^2

    discriminant = b^2 - a*c

    if discriminant ≥ 0
        d = √discriminant

        # t1 = (-b + d) / a
        t2 = (-b - d) / a

        if t2 > 0
            return t2
        end
    end

    return Inf   # no collision
end

"Flow ball for a time t"
function flow!(b::Particle, t, ::FreeFlow)
    b.x += b.v * t
end