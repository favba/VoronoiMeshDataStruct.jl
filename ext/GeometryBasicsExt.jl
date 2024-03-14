module GeometryBasicsExt

using VoronoiMeshDataStruct
using TensorsLite: Vec, norm
using TensorsLiteGeometry
using GeometryBasics: Polygon, Point2f, LineString, Line, TupleView

const PolType = Polygon{2, Float32, Point2f, LineString{2, Float32, Point2f, Base.ReinterpretArray{Line{2, Float32}, 1, Tuple{Point2f, Point2f}, TupleView{Tuple{Point2f, Point2f}, 2, 1, Vector{Point2f}}, false}}, Vector{LineString{2, Float32, Point2f, Base.ReinterpretArray{Line{2, Float32}, 1, Tuple{Point2f, Point2f}, TupleView{Tuple{Point2f, Point2f}, 2, 1, Vector{Point2f}}, false}}}}

function VoronoiMeshDataStruct.create_cells_polygons_periodic(vert_pos,cell_pos,verticesOnCell,x_period,y_period)

    cell_polygons = Vector{PolType}(undef,length(cell_pos))
    local_vertices = Vector{Point2f}(undef,VoronoiMeshDataStruct.max_length(eltype(verticesOnCell)))
    @inbounds for i in eachindex(cell_pos)
        cpos = cell_pos[i]
        l = 0
        for i_v in verticesOnCell[i]
            l+=1
            vpos = closest(cpos,vert_pos[i_v],x_period,y_period)
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

    vert_triangles = Vector{PolType}(undef,length(vert_pos))

    @inbounds for i in eachindex(vert_pos)
        vpos = vert_pos[i]
        ic1,ic2,ic3 = cellsOnVertex[i]

        cpos1 = closest(vpos,cell_pos[ic1],x_period,y_period)

        cpos2 = closest(vpos,cell_pos[ic2],x_period,y_period)

        cpos3 = closest(vpos,cell_pos[ic3],x_period,y_period)

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

    edge_quadrilaterals = Vector{PolType}(undef,length(edge_pos))

    @inbounds for i in eachindex(edge_pos)
        epos = edge_pos[i]
        iv1,iv2 = verticesOnEdge[i]
        ic1,ic2 = cellsOnEdge[i]

        cpos1 = closest(epos,cell_pos[ic1],x_period,y_period)

        vpos1 = closest(epos,vert_pos[iv1],x_period,y_period)

        cpos2 = closest(epos,cell_pos[ic2],x_period,y_period)

        vpos2 = closest(epos,vert_pos[iv2],x_period,y_period)

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
