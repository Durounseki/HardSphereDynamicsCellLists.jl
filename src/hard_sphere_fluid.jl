"""
Hard sphere fluid in N dimensions
"""
mutable struct HardSphereFluid{N,T}
	box::RectangularBox{N,T}
	particles::Vector{Particle{N,T}}
    cells::Dict{CartesianIndex{N}, Vector{Int64}}
    cell_L::T
end


#
# function HardSphereFluid(box::RectangularBox{N,T}, balls::Vector{MovableBall{N,T}}) where {N,T}
# 	partner1, partner2, collision_type, min_collision_time = find_collision(balls, box)
# 	current_time = zero(T)
#
# 	return HardSphereFluid{N,T}(box, balls, current_time, min_collision_time, partner1, partner2, collision_type)
#
# end



function HardSphereFluid{N,T}(box::RectangularBox{N,T}, num_particles, r, cell_L) where {N,T}

	particles = [Particle(zero(SVector{N,T}), zero(SVector{N,T}), r) for i in 1:num_particles]
    #Make a grid of the size of the maximum length
    maxL = findmax([findmax(box.lower)[1],findmax(box.upper)[1]])[1]
    cells = CellGrid(N, maxL, cell_L)

	return HardSphereFluid{N,T}(box, particles, cells, cell_L) #This line throws an error, instead of creating a variable of type HardsphereFluid it tries to run the function with missmatching types of arguments.
    #Let's try to move the cell_L property to box or declare the CellGrid as a struct

end



# HardSphereFluid(N, num_particles, r) = HardSphereFluid{N,Float64}(unit_hypercube(N, Float64), num_particles, r)
