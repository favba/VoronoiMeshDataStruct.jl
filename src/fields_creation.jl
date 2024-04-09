function compute_edge_normals_periodic(cpos::Vec2DxyArray,cellsOnEdge,xp::Number,yp::Number)
    n = similar(cpos,size(cellsOnEdge))

    @inbounds for i in eachindex(cellsOnEdge)
        ic1,ic2 = cellsOnEdge[i]
        cpos2 = cpos[ic2]
        cpos1 = closest(cpos2,cpos[ic1],xp,yp)
        n[i] = normalize(cpos2 - cpos1)
    end

    return n
end

function compute_edge_normals_on_sphere end

#This is only an approximation. Need to use great circles instead.
#function compute_edge_normals_on_sphere(cpos,cellsOnEdge)
#    n = similar(cpos,size(cellsOnEdge))
#
#    @inbounds for i in eachindex(cellsOnEdge)
#        ic1,ic2 = cellsOnEdge[i]
#        n[i] = normalize(cpos[ic2] - cpos[ic1])
#    end
#
#    return n
#end

compute_edge_normals(mesh::VoronoiMesh{true}) = compute_edge_normals_on_sphere(mesh.cells.position,mesh.edges.indices.cells)
compute_edge_normals(mesh::VoronoiMesh{false}) = compute_edge_normals_periodic(mesh.cells.position,mesh.edges.indices.cells,Float64(mesh.attributes[:x_period]),Float64(mesh.attributes[:y_period]))

function compute_edge_normals!(mesh::VoronoiMesh)
    mesh.edges.normalVectors = compute_edge_normals(mesh)
    return mesh.edges.normalVectors
end