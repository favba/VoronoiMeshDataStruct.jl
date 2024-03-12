module GeometryBasicsExt

using VoronoiMeshDataStruct
using TensorsLite: Vec, norm
using GeometryBasics: Polygon, Point2f, LineString, Line, TupleView

@inline function closest_point(p::Vec,points::Tuple{Vararg{<:Vec,N}}) where N
    @inline _,i = findmin(x->norm(x-p), points)
    r = @inbounds points[i]
    return r
end

@inline possible_positions(p::Vec,periods::Tuple{Vararg{<:Vec,N}}) where N = @inline map(+,ntuple(x->p,Val{N}()), periods)

const PolType = Polygon{2, Float32, Point2f, LineString{2, Float32, Point2f, Base.ReinterpretArray{Line{2, Float32}, 1, Tuple{Point2f, Point2f}, TupleView{Tuple{Point2f, Point2f}, 2, 1, Vector{Point2f}}, false}}, Vector{LineString{2, Float32, Point2f, Base.ReinterpretArray{Line{2, Float32}, 1, Tuple{Point2f, Point2f}, TupleView{Tuple{Point2f, Point2f}, 2, 1, Vector{Point2f}}, false}}}}

function VoronoiMeshDataStruct.create_cells_polygons_periodic(vert_pos,cell_pos,verticesOnCell,x_period,y_period)

    periods = (Vec(),
               Vec(x=x_period),Vec(x=-x_period),Vec(y=y_period),Vec(y=-y_period),
               Vec(x=x_period,y=y_period), Vec(x=x_period,y=-y_period), Vec(x=-x_period,y=y_period),Vec(x=-x_period,y=-y_period))

    cell_polygons = Vector{PolType}(undef,length(cell_pos))
    local_vertices = Vector{Point2f}(undef,VoronoiMeshDataStruct.max_length(eltype(verticesOnCell)))
    @inbounds for i in eachindex(cell_pos)
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

function VoronoiMeshDataStruct.create_dual_triangles_periodic(vert_pos,cell_pos,cellsOnVertex,x_period,y_period)
    periods = (Vec(),
               Vec(x=x_period),Vec(x=-x_period),Vec(y=y_period),Vec(y=-y_period),
               Vec(x=x_period,y=y_period), Vec(x=x_period,y=-y_period), Vec(x=-x_period,y=y_period),Vec(x=-x_period,y=-y_period))

    vert_triangles = Vector{PolType}(undef,length(vert_pos))

    @inbounds for i in eachindex(vert_pos)
        vpos = vert_pos[i]
        ic1,ic2,ic3 = cellsOnVertex[i]

        possible_cpos = possible_positions(cell_pos[ic1], periods) 
        cpos1 = closest_point(vpos,possible_cpos)

        possible_cpos = possible_positions(cell_pos[ic2], periods) 
        cpos2 = closest_point(vpos,possible_cpos)

        possible_cpos = possible_positions(cell_pos[ic3], periods) 
        cpos3 = closest_point(vpos,possible_cpos)

        vert_triangles[i] = Polygon([Point2f(cpos1.x,cpos1.y),Point2f(cpos2.x,cpos2.y),Point2f(cpos3.x,cpos3.y)])
    end

    return vert_triangles
end

function VoronoiMeshDataStruct.create_dual_triangles(mesh::VoronoiMesh)
    if mesh.attributes[:is_periodic]::String == "YES"
        return create_dual_triangles_periodic(mesh.vertices.position, mesh.cells.position, mesh.vertices.indices.cells, mesh.attributes[:x_period], mesh.attributes[:y_period]) 
    else
        error("Not yet implemented for non-periodic meshes")
    end
end

function VoronoiMeshDataStruct.create_edge_quadrilaterals_periodic(edge_pos,vert_pos,cell_pos,verticesOnEdge,cellsOnEdge,x_period,y_period)
    periods = (Vec(),
               Vec(x=x_period),Vec(x=-x_period),Vec(y=y_period),Vec(y=-y_period),
               Vec(x=x_period,y=y_period), Vec(x=x_period,y=-y_period), Vec(x=-x_period,y=y_period),Vec(x=-x_period,y=-y_period))

    edge_quadrilaterals = Vector{PolType}(undef,length(edge_pos))

    @inbounds for i in eachindex(edge_pos)
        epos = edge_pos[i]
        iv1,iv2 = verticesOnEdge[i]
        ic1,ic2 = cellsOnEdge[i]

        possible_pos = possible_positions(cell_pos[ic1], periods) 
        cpos1 = closest_point(epos,possible_pos)

        possible_pos = possible_positions(vert_pos[iv1], periods) 
        vpos1 = closest_point(epos,possible_pos)

        possible_pos = possible_positions(cell_pos[ic2], periods) 
        cpos2 = closest_point(epos,possible_pos)

        possible_pos = possible_positions(vert_pos[iv2], periods) 
        vpos2 = closest_point(epos,possible_pos)


        edge_quadrilaterals[i] = Polygon([Point2f(cpos1.x,cpos1.y),Point2f(vpos1.x,vpos1.y),Point2f(cpos2.x,cpos2.y),Point2f(vpos2.x,vpos2.y)])
    end

    return edge_quadrilaterals
end

function VoronoiMeshDataStruct.create_edge_quadrilaterals(mesh::VoronoiMesh)
    if mesh.attributes[:is_periodic]::String == "YES"
        return create_edge_quadrilaterals_periodic(mesh.edges.position, mesh.vertices.position, mesh.cells.position, mesh.edges.indices.vertices, mesh.edges.indices.cells, mesh.attributes[:x_period], mesh.attributes[:y_period]) 
    else
        error("Not yet implemented for non-periodic meshes")
    end
end

end
