"""
	Wall{N,T}

`N`-dimensional hyperplane.

- `p`: point on the plane
- `n`: normal vector
- `Temp`: temperature
"""
mutable struct Wall{N,T}
    p::SVector{N,T}
    n::SVector{N,T}
    Temp::T
end

function normal(p::Wall)
    p.n
end


"""
"Rectangular" box in N dimensions.

`lower` and `upper` contain the lower and upper bounds in each dimension.
"""
mutable struct RectangularBox{N,T}
    lower::SVector{N,T}
    upper::SVector{N,T}
    walls::Vector{Wall{N,T}}
    Temp::T
    cells::Dict{CartesianIndex{N}, Vector{Int64}}
    cell_L::T

    function RectangularBox{N,T}(lower::SVector{N,T}, upper::SVector{N,T},Temp::T,cell_L::T) where {N,T}

        z = zero(SVector{N,T})

        walls = Wall{N,T}[]

        for i in 1:N
            n = setindex(z, -1, i)  # unit vector with one non-zero component

            push!(walls, Wall(lower, n, Temp))
            push!(walls, Wall(upper, -n, Temp))
        end

        Lmax = findmax([findmax(lower)[1],findmax(upper)[1]])[1]
        Lmin = findmin([findmin(lower)[1],findmin(upper)[1]])[1]
        
        cells = CellGrid(N, Lmin, Lmax, cell_L)

        new{N,T}(lower, upper, walls, Temp, cells, cell_L)
    end
end

RectangularBox(lower::SVector{N,T}, upper::SVector{N,T}, Temp::T,cell_L::T) where {N,T} = RectangularBox{N,T}(lower, upper, Temp, cell_L)

# function unit_hypercube(N, T)
# 	return RectangularBox(-0.5 .* ones(SVector{N,T}), +0.5 .* ones(SVector{N,T}), Temp)
# end