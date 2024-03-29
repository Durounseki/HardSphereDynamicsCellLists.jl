@inline function neighbors(d::Int)
    n = CartesianIndices((fill(-1:1, d)...,))
    return n[1:fld(length(n), 2)]
end

@inline function CellGrid(d::Int, Lmin, Lmax, r₀)
    indices = CartesianIndices((fill((Int(floor(Lmin/r₀))):Int(floor(Lmax/r₀)), d)...,))
    grid = Dict{CartesianIndex{d}, Vector{Int}}()
    for ci in indices
        grid[ci]=[]
    end
    return grid
end

"""

Generate list of particles in each cell
`positions`: array of particles positions
`r₀`: size of cell
`cells`: Dictionary of empty arrays with cartesian indices as keys

Add the indices of particles contained in each cell to `cells[(i₁,i₂,...,iₙ)`]

"""

function CellList(positions, box::RectangularBox{N,T}) where {N,T}
    n=length(positions)
    for i in 1:n
        push!(box.cells[CartesianIndex(Int.(floor.(positions[i] ./ box.cell_L))...)],i)
    end
end

function ghostCellsTB(Lx,Ly,r₀)
    tbBoundary=[CartesianIndex(i,Int(floor(Ly/r₀))-1) for i in 0:Int(floor(Lx/r₀))-1]
    append!(tbBoundary,[CartesianIndex(i,0) for i in 0:Int(floor(Lx/r₀))-1])
    return tbBoundary
end

function ghostCellsRL(Lx,Ly,r₀)
    rlBoundary=[CartesianIndex(Int(floor(Lx/r₀))-1,i) for i in 0:Int(floor(Ly/r₀))-1]
    append!(rlBoundary,[CartesianIndex(0,i) for i in 0:Int(floor(Ly/r₀))-1])
    return rlBoundary
end


@inline function cellNeighbors!(ps,ds, is, p, r₀)
    for (k, i) in enumerate(is[1:(end-1)])
        for j in is[(k+1):end]
            r=norm(p[i]-p[j])
            if r <= r₀
                # println(i,j,p[i],p[j],norm(p[i]-p[j]))
                ps[i][j]=1
                ds[i][j]=r
                ps[j][i]=1
                ds[j][i]=r
            else
                ps[i][j]=0
                ds[i][j]=0.0
                ps[j][i]=0
                ds[j][i]=0.0
            end
        end
        # ps[i][i]=(0,0.0)
    end
end

@inline function cellNeighbors!(ps, ds, is, js, p, r₀)
    for i in is
        for j in js
            r = norm(p[i]-p[j])
            if r <= r₀
                # println(i,j,p[i],p[j],norm(p[i]-p[j]))
                ps[i][j]=1
                ds[i][j]=r
                ps[j][i]=1
                ds[j][i]=r
            else
                ps[i][j]=0
                ds[i][j]=0.0
                ps[j][i]=0
                ds[j][i]=0.0
            end
        end
        # ps[i][i]=(0,0.0)
    end
end

"""

Calculate the nearest neighbour on one cell

`particles`: collection of particles
`indices`: denotes the particles contained in one cell

returns the smallest distance between pairs of particles and their indices

"""
@inline function nearestNeighbor!(particles,indices)
    r_min = Inf
    partner1 = -1
    partner2 = -1
    for (k, i) in enumerate(indices[1:(end-1)])
        for j in indices[(k+1):end]
            r=norm(particles[i].x-particles[j].x)
            if r < r_min
                r_min = r
                partner1 = i
                partner2 = j
            end
        end
    end
    return r_min, partner1, partner2
end

##Collision time register

# function pairCollisionTimes(N)

#     collision_times = Dict{NTuple{N,Int},Float64}()
#     for i in 1:N
#         for j in 

# end


#Next colliding pair within one cell
"""

Calculate the nearest neighbour on one cell

`particles`: collection of particles
`indices`: denotes the particles contained in one cell

returns the smallest distance between pairs of particles and their indices

"""
@inline function collidingPair!(particles::Vector{Particle{N,T}}, flow, indices) where{N,T}
    
    p_p_collision_time = 0.0
    
    for (k, i) in enumerate(indices[1:(end-1)])
        for j in indices[(k+1):end]
            p_p_collision_time = collision_time(particles[i], particles[j], flow)
            if p_p_collision_time < particles[i].col_time
                particles[i].col_time = p_p_collision_time
                particles[i].col_pair = j
                # if particles[i].col_type != :disc_collision
                #     particles[i].col_type = :disc_collision
                # end
            end
            if p_p_collision_time < particles[j].col_time
                particles[j].col_time = p_p_collision_time
                particles[j].col_pair = i
                # if particles[i].col_type != :disc_collision
                #     particles[i].col_type = :disc_collision
                # end
            end
        end
    end

end

#####

#Check for next colliding pair in neighboring cells

"""

If there are no more than one particle in one cell, calculate the nearest neighbour in the adjacent cells

`particles`: collection of particles
`center_index`: denotes the particle contained within the center cell
`indices`: denotes the particles contained in the adjacent cells

returns the smallest distance between pairs of particles and their indices

"""
@inline function collidingPair!(particles, flow, center_indices, indices)
    
    p_p_collision_time = 0.0
    
    for i in center_indices
        for j in indices
            p_p_collision_time = collision_time(particles[i], particles[j], flow)
            if p_p_collision_time < particles[i].col_time
                particles[i].col_time = p_p_collision_time
                particles[i].col_pair = j
                # if particles[i].col_type != :disc_collision
                #     particles[i].col_type = :disc_collision
                # end
            end
            if p_p_collision_time < particles[j].col_time
                particles[j].col_time = p_p_collision_time
                particles[j].col_pair = i
                # if particles[i].col_type != :disc_collision
                #     particles[i].col_type = :disc_collision
                # end
            end
        end
    end
end

######

function near_neighbors(cells, P, r₀, FRN, ds)
    offsets = neighbors(length(P[1]))
    # Iterate over non-empty cells
    for (cell, is) in cells
        # Pairs of points within the cell
        cellNeighbors!(FRN, ds, is, P, r₀)
        # Pairs of points with non-empty neighboring cells
        for offset in offsets
            neigh_cell = cell + offset
            if haskey(cells, neigh_cell)
                @inbounds js = cells[neigh_cell]
                cellNeighbors!(FRN, ds, is, js, P, r₀)
            end
        end
    end
end


# #Calculate next colliding pair

# """

# Calculate the pair of particles (i,j) with the smallest interparticle distance rᵢⱼ

# `particles`: collection of particles
# `cells`: dictionary of indices with cartesian indices as keys denoting the cells in the grid

# Returns the indices of the next colliding pair

# """
# function nextCollidingPair(cells, particles, flow, min_collision_times)
    
#     # min_collision_time = Inf
#     # partner1 = -1
#     # partner2 = -1

#     offsets = neighbors(length(particles[1].x))

#     for (cell, center_indices) in cells
#         # Iterate over non-empty cells
#         if length(center_indices) > 0

#             # Check for next colliding pair within the same cell

#             collidingPair!(particles, flow, center_indices, min_collision_times)

#             # Check for next colliding pair in neighboring cells

#             for offset in offsets
#                 neigh_cell = cell + offset
#                 if haskey(cells, neigh_cell)
#                     @inbounds indices = cells[neigh_cell]
#                     collidingPair!(particles, flow, center_indices, indices, min_collision_times)
#                 end
#             end
#         end
#     end

#     time_to_collision, pair = findmin(min_collision_times)

#     #If all the cells contain a single particle, the system is probably dilute, or the cells are very small.
#     #In this case we could increase the size of the cell or we can try to iterate over all the pairs.
#     # if findmax([length(i) for (c,i) in cells])[1] < 2

#     #     for i in 1:length(particles)
#     #         for j in i+1:length(particles)
                
#     #             t = collision_time(particles[i], particles[j], flow)

#     #             if t < min_collision_time
#     #                 partner1 = i
#     #                 partner2 = j
#     #                 min_collision_time = t
#     #             end

#     #         end
#     #     end

#     # else
#     #     for (cell, center_indices) in cells
#     #         # Iterate over non-empty cells
#     #         if length(center_indices) > 0

#     #             # Check for next colliding pair within the same cell

#     #             collidingPair!(particles, flow, center_indices, min_collision_time, partner1, partner2)

#     #             # Check for next colliding pair in neighboring cells

#     #             for offset in offsets
#     #                 neigh_cell = cell + offset
#     #                 if haskey(cells, neigh_cell)
#     #                     @inbounds indices = cells[neigh_cell]
#     #                     collidingPair!(particles, flow, center_indices, indices, min_collision_time, partner1, partner2)
#     #                 end
#     #             end
#     #         end
#     #     end
#     # end

#     return time_to_collision, pair[1], pair[2]

# end


# ######


# """

# Calculate the pair of particles (i,j) with the smallest interparticle distance rᵢⱼ

# `particles`: collection of particles
# `cells`: dictionary of indices with cartesian indices as keys denoting the cells in the grid

# Returns the indices of the next colliding pair

# """
# function nextCollidingPairOld(cells, particles)
    
#     r_min = Inf
#     partner1 = -1
#     partner2 = -1

#     offsets = neighbors(length(particles[1].x))

#     #If all the cells contain a single particle, the system is probably dilute, or the cells are very small.
#     #In this case we could increase the size of the cell or we can try to iterate over all the pairs.
#     if findmax([length(i) for (c,i) in cells])[1] < 2

#         for i in 1:length(particles)
#             for j in i+1:length(particles)
#                 r = norm(particles[i].x-particles[j].x)
#                 if r < r_min
#                     r_min = r
#                     partner1 = i
#                     partner2 = j
#                 end
#             end
#         end

#     else
#         for (cell, indices) in cells
#             # Iterate over non-empty cells
#             if length(indices) > 0

#                 # Pairs of points within the cell
#                 if length(indices) > 1
#                     r, i, j = nearestNeighbor!(particles,indices)
#                 else
#                     r = Inf
#                     i = -1
#                     j = -1
#                     # If cell contains less than two particles, check for pairs of points in neighboring cells
#                     for offset in offsets
#                         neigh_cell = cell + offset
#                         if haskey(cells, neigh_cell)
#                             @inbounds indices2 = cells[neigh_cell]
#                             #center_index = indices[1]
#                             r_temp, i_temp, j_temp = nearestNeighbor!(particles, indices[1], indices2)
#                             if r_temp < r
#                                 r = r_temp
#                                 i = i_temp
#                                 j = j_temp
#                             end
#                         end
#                     end
#                 end

#                 if r < r_min
#                     r_min = r
#                     partner1 = i
#                     partner2 = j
#                 end
                
#             end
#         end
#     end

#     return partner1, partner2
#     # #Next colliding pair
#     # particles[i].nn = 1
#     # particles[j].nn = 1

# end

function emptyCells(cells)
    for (ci,cell) in cells
        while length(cell)>0
            pop!(cell)
        end
    end
end

function nextCollidingPair(cells, particles, flow)

    offsets = neighbors(length(particles[1].x))
    neigh_cell = CartesianIndex(Int.(zero(particles[1].x))...)

    for (cell, indices) in cells
        collidingPair!(particles, flow, indices)
        for offset in offsets
            neigh_cell = cell + offset
            if haskey(cells, neigh_cell)
                @inbounds indices2 = cells[neigh_cell]
                collidingPair!(particles, flow, indices, indices2)
            end
        end
    end

end