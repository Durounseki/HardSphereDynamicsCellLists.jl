"""
Hard sphere fluid in N dimensions
"""
mutable struct HardSphereFluid{N,T}
	box::RectangularBox{N,T}
	particles::Vector{Particle{N,T}}
    # cells::Dict{CartesianIndex{N}, Vector{Int64}}
end


#
# function HardSphereFluid(box::RectangularBox{N,T}, balls::Vector{MovableBall{N,T}}) where {N,T}
# 	partner1, partner2, collision_type, min_collision_time = find_collision(balls, box)
# 	current_time = zero(T)
#
# 	return HardSphereFluid{N,T}(box, balls, current_time, min_collision_time, partner1, partner2, collision_type)
#
# end



function HardSphereFluid{N,T}(box::RectangularBox{N,T}, num_particles, r) where {N,T} #We still need a better solution so that we can use box

	particles = [Particle(zero(SVector{N,T}), zero(SVector{N,T}), r) for i in 1:num_particles]
    #Make a grid of the size of the maximum length
    # maxL = findmax([findmax(box.lower)[1],findmax(box.upper)[1]])[1]
    # cells = CellGrid(N, maxL, cell_L)

	return HardSphereFluid{N,T}(box, particles)

end



# HardSphereFluid(N, num_particles, r) = HardSphereFluid{N,Float64}(unit_hypercube(N, Float64), num_particles, r)
