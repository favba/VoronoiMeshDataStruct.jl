function graph_partition(cellsOnCell::Vector{<:ImmutableVector},nEdges::Integer)
    nCells = length(cellsOnCell)
    outIO = IOBuffer()
    println(outIO, nCells,' ',nEdges)

    @inbounds for i in eachindex(cellsOnCell)
        map(x->print(outIO,x,' '),cellsOnCell[i])
        skip(outIO,-1)
        println(outIO)
    end
    seekstart(outIO)
    return outIO
end

graph_partition(cells::Union{CellBase,CellInfo},edges::Union{EdgeBase,EdgeInfo}) = graph_partition(cells.indices.cells,edges.n)

graph_partition(mesh::VoronoiMesh) = graph_partition(mesh.cells,mesh.edges)

function find_obtuse_triangles_periodic(cpos,cellsOnVertex,xp::Number,yp::Number)
    r = Int[]
    lk = ReentrantLock()
    @parallel for v in eachindex(cellsOnVertex)
        @inbounds begin
        c1,c2,c3 = cellsOnVertex[v]
        c1pos = cpos[c1]
        c2pos = closest(c1pos,cpos[c2],xp,yp)
        c3pos = closest(c1pos,cpos[c3],xp,yp)
        if is_obtuse(c1pos,c2pos,c3pos)
            lock(lk) do
                push!(r,v)
            end
        end
        end #inbounds
    end
    return r
end

find_obtuse_triangles(vertices::VertexBase{false},cells::CellBase{false},xp::Number,yp::Number) = find_obtuse_triangles_periodic(cells.position,vertices.indices.cells,xp,yp)

find_obtuse_triangles(mesh::VoronoiMesh{false}) = find_obtuse_triangles(mesh.vertices.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

"""
    periodic_edges_mask(dcEdge::Vector,cellsOnEdge,cell_positions)
    periodic_edges_mask(mesh::VoronoiMesh)

Returns a BitArray that masks periodic edges (interior edges = true, boundary edges = false).
"""
function periodic_edges_mask(dc,cellsOnEdge,c_pos)
    mask = BitArray(undef,length(dc))
    @inbounds for e in eachindex(cellsOnEdge)
        dc_e = dc[e]
        c1,c2 = cellsOnEdge[e]
        dc_2 = norm(c_pos[c2] - c_pos[c1])
        mask[e] = dc_e ≈ dc_2
    end
    return mask
end

periodic_edges_mask(mesh::VoronoiMesh) = periodic_edges_mask(mesh.edges.dc, mesh.edges.indices.cells, mesh.cells.position)