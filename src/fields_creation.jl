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

#Since the edges are always equidistant to both cell positions, the vector connecting them will be parallel to the arc tangent at the edge, right??
function compute_edge_normals_on_sphere(cpos,cellsOnEdge)
    n = similar(cpos,size(cellsOnEdge))

    @inbounds for i in eachindex(cellsOnEdge)
        ic1,ic2 = cellsOnEdge[i]
        n[i] = normalize(cpos[ic2] - cpos[ic1])
    end

    return n
end

compute_edge_normals(mesh::VoronoiMesh{true}) = compute_edge_normals_on_sphere(mesh.cells.position,mesh.edges.indices.cells)
compute_edge_normals(mesh::VoronoiMesh{false}) = compute_edge_normals_periodic(mesh.cells.position,mesh.edges.indices.cells,Float64(mesh.attributes[:x_period]),Float64(mesh.attributes[:y_period]))

function compute_edge_normals!(mesh::VoronoiMesh)
    mesh.edges.normalVectors = compute_edge_normals(mesh)
    return mesh.edges.normalVectors
end

function compute_area_triangles_periodic!(output,cpos,cellsOnVertex,xp::Number,yp::Number)
    @inbounds for v in eachindex(cellsOnVertex)
        c1,c2,c3 = cellsOnVertex[v]
        c1_pos = cpos[c1]
        c2_pos = closest(c1_pos,cpos[c2],xp,yp)
        c3_pos = closest(c1_pos,cpos[c3],xp,yp)
        output[v] = area(c1_pos,c2_pos,c3_pos)
    end
    return output
end

function compute_area_triangles_periodic(cpos,cellsOnVertex,xp,yp)
    output = Vector{nonzero_eltype(eltype(cpos))}(undef,length(cellsOnVertex))
    return compute_area_triangles_periodic!(output,cpos,cellsOnVertex,xp,yp)
end

compute_area_triangles(mesh::VoronoiMesh{false}) = compute_area_triangles_periodic(mesh.cells.position,mesh.vertices.indices.cells,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_area_triangles!(output::AbstractVector,mesh::VoronoiMesh{false})
    length(output) == mesh.vertices.n || throw(DomainError("Output array doesn't have correct length.\nArray length: "*length(output)*".\nNumber of vertices:"*mesh.vertices.n))
    return compute_area_triangles_periodic!(output,mesh.cells.position,mesh.vertices.indices.cells,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)
end

compute_area_triangles!(mesh::VoronoiMesh{false}) = compute_area_triangles!(mesh.vertices.area,mesh)

function compute_kite_areas_periodic!(output,cpos,vpos,cellsOnVertex,xp,yp)
  @inbounds for v in eachindex(cellsOnVertex)
        v_pos = vpos[v]
        c1,c2,c3 = cellsOnVertex[v]
        c1_pos = closest(v_pos,cpos[c1],xp,yp)
        c2_pos = closest(v_pos,cpos[c2],xp,yp)
        c3_pos = closest(v_pos,cpos[c3],xp,yp)
        
        c12_pos = 0.5*(c1_pos + c2_pos)
        c23_pos = 0.5*(c2_pos + c3_pos)
        c31_pos = 0.5*(c3_pos + c1_pos)

        a1 = area(c1_pos,c12_pos,v_pos,c31_pos)
        a2 = area(c12_pos,c2_pos,c23_pos,v_pos)
        a3 = area(c23_pos,c3_pos,c31_pos,v_pos)
        output[v] = (a1,a2,a3)
    end
    return output
end

function compute_kite_areas_periodic(cpos,vpos,cellsOnVertex,xp,yp)
    output = Vector{NTuple{3,nonzero_eltype(eltype(cpos))}}(undef,length(cellsOnVertex))
    return compute_kite_areas_periodic!(output,cpos,vpos,cellsOnVertex,xp,yp)
end

compute_kite_areas(mesh::VoronoiMesh{false}) = compute_kite_areas_periodic(mesh.cells.position,mesh.vertices.position,mesh.vertices.indices.cells,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_kite_areas!(output::AbstractVector,mesh::VoronoiMesh{false})
    length(output) == mesh.vertices.n || throw(DomainError("Output array doesn't have correct length.\nArray length: "*length(output)*".\nNumber of vertices:"*mesh.vertices.n))
    return compute_kite_areas_periodic!(output,mesh.cells.position,mesh.vertices.position,mesh.vertices.indices.cells,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)
end
