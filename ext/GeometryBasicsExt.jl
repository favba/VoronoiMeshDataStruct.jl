module GeometryBasicsExt

using VoronoiMeshDataStruct
using TensorsLite: Vec, norm
using GeometryBasics: Polygon, Point2f

@inline function closest_point(p::Vec,points::Tuple{Vararg{<:Vec,N}}) where N
    @inline _,i = findmin(x->norm(x-p), points)
    r = @inbounds points[i]
    return r
end

@inline possible_positions(p::Vec,periods::Tuple{Vararg{<:Vec,N}}) where N = @inline map(+,ntuple(x->p,Val{N}()), periods)

function VoronoiMeshDataStruct.create_cells_polygons_periodic(vert_pos,cell_pos,verticesOnCell,x_period,y_period)

    periods = (Vec(),
               Vec(x=x_period),Vec(x=-x_period),Vec(y=y_period),Vec(y=-y_period),
               Vec(x=x_period,y=y_period), Vec(x=x_period,y=-y_period), Vec(x=-x_period,y=y_period),Vec(x=-x_period,y=-y_period))

    local_vertices = Vector{Point2f}(undef,VoronoiMeshDataStruct.max_length(eltype(verticesOnCell)))
    l = 0
    @inbounds for i_v in verticesOnCell[1]
        l+=1
        possible_vpos = possible_positions(vert_pos[i_v], periods) 
        cpos = cell_pos[1]
        vpos = closest_point(cpos,possible_vpos)
        local_vertices[l] = Point2f(vpos.x,vpos.y)
    end

    polygon1 = Polygon(local_vertices[1:l])
    cell_polygons = Vector{typeof(polygon1)}(undef,length(cell_pos))
    cell_polygons[1] = polygon1

    @inbounds for i in eachindex(cell_pos)[2:end]
        cpos = cell_pos[i]
        l = 0
        for i_v in verticesOnCell[i]
            l+=1
            possible_vpos = possible_positions(vert_pos[i_v], periods) 
            vpos = closest_point(cpos,possible_vpos)
            local_vertices[l] = Point2f(vpos.x,vpos.y)
        end
        cell_polygons[i] = Polygon(local_vertices[1:l])
    end

    return cell_polygons
end

function VoronoiMeshDataStruct.create_cells_polygons(mesh::VoronoiMesh)
    if mesh.attributes[:is_periodic]::String == "YES"
        return create_cells_polygons_periodic(mesh.vertices.position, mesh.cells.position, mesh.cells.indices.vertices, mesh.attributes[:x_period], mesh.attributes[:y_period]) 
    else
        error("Not yet implemented for non-periodic meshes")
    end
end

end
