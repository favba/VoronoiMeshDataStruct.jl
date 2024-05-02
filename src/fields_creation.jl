function compute_edge_normals_periodic!(n::Vec2DxyArray,cpos::Vec2DxyArray,cellsOnEdge,xp::Number,yp::Number)

    @inbounds for i in eachindex(cellsOnEdge)
        ic1,ic2 = cellsOnEdge[i]
        cpos2 = cpos[ic2]
        cpos1 = closest(cpos2,cpos[ic1],xp,yp)
        n[i] = normalize(cpos2 - cpos1)
    end

    return n
end

function compute_edge_normals_periodic(cpos::Vec2DxyArray,cellsOnEdge,xp::Number,yp::Number)
    n = similar(cpos,size(cellsOnEdge))
    return compute_edge_normals_periodic!(n,cpos,cellsOnEdge,xp,yp)
end

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

function compute_edge_normals!(output,mesh::VoronoiMesh{false})
    length(output) == mesh.edges.n || throw(DomainError("Output array doesn't have correct length.\nArray length: "*length(output)*".\nNumber of edges:"*mesh.edges.n))
    compute_edge_normals_periodic!(output,mesh.cells.position,mesh.edges.indices.cells,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)
end

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
        
        #Those should be the edges positions
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

function compute_area_cells_periodic!(output,vpos,verticesOnCell,xp::Number,yp::Number)
    @inbounds for c in eachindex(verticesOnCell)
        output[c] = area(vpos,verticesOnCell[c],xp,yp)
    end
    return output
end

function compute_area_cells_periodic(vpos,verticesOnCell,xp::Number,yp::Number)
    output = Vector{nonzero_eltype(eltype(vpos))}(undef,length(verticesOnCell))
    return compute_area_cells_periodic!(output,vpos,verticesOnCell,xp,yp)
end

compute_area_cells(mesh::VoronoiMesh{false}) = compute_area_cells_periodic(mesh.vertices.position,mesh.cells.indices.vertices,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_area_cells!(output::AbstractVector,mesh::VoronoiMesh{false})
    length(output) == mesh.cells.n || throw(DomainError("Output array doesn't have correct length.\nArray length: "*length(output)*".\nNumber of cells:"*mesh.cells.n))
    return compute_area_cells_periodic!(output,mesh.vertices.position,mesh.cells.indices.vertices,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)
end

compute_area_cells!(mesh::VoronoiMesh{false}) = compute_area_cells!(mesh.cells.area,mesh)

function compute_dcEdge_periodic!(output,cpos,cellsOnEdge,xp::Number,yp::Number)
    @inbounds for e in eachindex(cellsOnEdge)
        c1,c2 = cellsOnEdge[e]
        c1pos = cpos[c1]
        c2pos = closest(c1pos,cpos[c2],xp,yp)
        output[e] = norm(c2pos-c1pos)
    end
    return output
end

function compute_dcEdge!(output,mesh::VoronoiMesh{false}) 
    length(output) == mesh.edges.n || throw(DomainError("Output array doesn't have correct length.\nArray length: "*length(output)*".\nNumber of edges:"*mesh.edges.n))
    return compute_dcEdge_periodic!(output,mesh.cells.position,mesh.edges.indices.cells,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)
end

function compute_dcEdge_periodic(cpos,cellsOnEdge,xp::Number,yp::Number)
    output = Vector{nonzero_eltype(eltype(cpos))}(undef,length(cellsOnEdge))
    return compute_dcEdge_periodic!(output,cpos,cellsOnEdge,xp,yp)
end

compute_dcEdge(mesh::VoronoiMesh{false}) = compute_dcEdge_periodic(mesh.cells.position,mesh.edges.indices.cells,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_dcEdge!(mesh::VoronoiMesh) 
    if isdefined(mesh.edges,:dc)
        return compute_dcEdge!(mesh.edges.dc,mesh)
    else
        mesh.edges.dc = compute_dcEdge(mesh)
        return mesh.edges.dc
    end
end

function compute_dvEdge_periodic!(output,vpos,verticesOnEdge,xp::Number,yp::Number)
    @inbounds for e in eachindex(verticesOnEdge)
        v1,v2 = verticesOnEdge[e]
        v1pos = vpos[v1]
        v2pos = closest(v1pos,vpos[v2],xp,yp)
        output[e] = norm(v2pos-v1pos)
    end
    return output
end

function compute_dvEdge!(output,mesh::VoronoiMesh{false}) 
    length(output) == mesh.edges.n || throw(DomainError("Output array doesn't have correct length.\nArray length: "*length(output)*".\nNumber of edges:"*mesh.edges.n))
    return compute_dvEdge_periodic!(output,mesh.vertices.position,mesh.edges.indices.vertices,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)
end

function compute_dvEdge_periodic(vpos,verticesOnEdge,xp::Number,yp::Number)
    output = Vector{nonzero_eltype(eltype(vpos))}(undef,length(verticesOnEdge))
    return compute_dvEdge_periodic!(output,vpos,verticesOnEdge,xp,yp)
end

compute_dvEdge(mesh::VoronoiMesh{false}) = compute_dvEdge_periodic(mesh.vertices.position,mesh.edges.indices.vertices,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_dvEdge!(mesh::VoronoiMesh) 
    if isdefined(mesh.edges,:dv)
        return compute_dvEdge!(mesh.edges.dv,mesh)
    else
        mesh.edges.dv = compute_dvEdge(mesh)
        return mesh.edges.dv
    end
end

function compute_angleEdge_periodic!(output,cpos,cellsOnEdge,xp::Number,yp::Number)
  @inbounds for i in eachindex(cellsOnEdge)
        ic1,ic2 = cellsOnEdge[i]
        cpos2 = cpos[ic2]
        cpos1 = closest(cpos2,cpos[ic1],xp,yp)
        output[i] = acos(normalize(cpos2 - cpos1) ⋅ 𝐢)
    end
    return output
end

function compute_angleEdge!(output,mesh::VoronoiMesh{false}) 
    length(output) == mesh.edges.n || throw(DomainError("Output array doesn't have correct length.\nArray length: "*length(output)*".\nNumber of edges:"*mesh.edges.n))
    return compute_angleEdge_periodic!(output,mesh.cells.position,mesh.edges.indices.cells,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)
end

function compute_angleEdge_periodic(cpos,cellsOnEdge,xp::Number,yp::Number)
    output = Vector{nonzero_eltype(eltype(cpos))}(undef,length(cellsOnEdge))
    return compute_angleEdge_periodic!(output,cpos,cellsOnEdge,xp,yp)
end

compute_angleEdge(mesh::VoronoiMesh{false}) = compute_angleEdge_periodic(mesh.cells.position,mesh.edges.indices.cells,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_angleEdge!(mesh::VoronoiMesh) 
    if isdefined(mesh.edges,:angle)
        return compute_angleEdge!(mesh.edges.angle,mesh)
    else
        mesh.edges.dc = compute_angleEdge(mesh)
        return mesh.edges.angle
    end
end

function select_shared_vertex(next_vs, previous_vs)
    v = zero(eltype(next_vs))
    nv1,nv2 = next_vs
    pv1,pv2 = previous_vs

    if ((nv1 == pv1) || (nv1 == pv2))
        v = nv1
    elseif ((nv2 == pv1) || (nv2 == pv2))
        v = nv2
    else
        error("No common vertex found")
    end
    return v
end

function select_kite_area(kiteAreaOnVertex,cellsOnVertex,v,c)
    areas = kiteAreaOnVertex[v]
    cells = cellsOnVertex[v]
    if cells[1] == c
        return areas[1]
    elseif cells[2] == c
        return areas[2]
    elseif cells[3] == c
        return areas[3]
    else
        error("Vertex $v doesn't belong to cell $c")
    end
end

function sign_edge((c1,c2),c)
    if c == c1
        return 1
    elseif c == c2
        return -1
    else
        error("Edge doesn't belong to cell $c")
    end
end

function compute_weightsOnEdge_trisk!(weightsOnEdge,verticesOnEdge,cellsOnEdge,edgesOnEdge,dcEdge,dvEdge,kiteAreasOnVertex,cellsOnVertex,nEdgesOnCell,areaCell)
    for e in eachindex(edgesOnEdge)
        c1,c2 = cellsOnEdge[e]
        inds_e = edgesOnEdge[e]
        inv_de = inv(dcEdge[e])

        nEdges_c1 = nEdgesOnCell[c1]
        inv_a_c1 = inv(areaCell[c1])
        previous_vs = verticesOnEdge[e]
        R = zero(eltype(weightsOnEdge))
        for i in 1:(nEdges_c1-1)
            next_e = inds_e[i]
            next_vs = verticesOnEdge[next_e]
            v = select_shared_vertex(next_vs, previous_vs)
            Avi = select_kite_area(kiteAreasOnVertex,cellsOnVertex,v,c1)
            R += inv_a_c1*Avi
            weightsOnEdge[i,e] = sign_edge(cellsOnEdge[next_e],c1)*inv_de*dvEdge[next_e]*(0.5 - R)
            previous_vs = next_vs
        end

        nEdges_c2 = nEdgesOnCell[c2]
        inv_a_c2 = inv(areaCell[c2])
        R = zero(eltype(weightsOnEdge))
        for i in nEdges_c1:(nEdges_c1 + nEdges_c2 - 2)
            next_e = inds_e[i]
            next_vs = verticesOnEdge[next_e]
            v = select_shared_vertex(next_vs, previous_vs)
            Avi = select_kite_area(kiteAreasOnVertex,cellsOnVertex,v,c2)
            R += inv_a_c2*Avi
            weightsOnEdge[i,e] = (-sign_edge(cellsOnEdge[next_e],c2))*inv_de*dvEdge[next_e]*(0.5 - R)
            previous_vs = next_vs
        end
    end
    return weightsOnEdge
end

function compute_weightsOnEdge_trisk(verticesOnEdge,cellsOnEdge,edgesOnEdge,dcEdge,dvEdge,kiteAreasOnVertex,cellsOnVertex,nEdgesOnCell,areaCell)
    T = eltype(eltype(kiteAreasOnVertex))
    nedges_max = max_length(eltype(edgesOnEdge))
    nedges = length(verticesOnEdge)
    weightsOnEdge = zeros(T,nedges_max,nedges)
    return compute_weightsOnEdge_trisk!(weightsOnEdge,verticesOnEdge,cellsOnEdge,edgesOnEdge,dcEdge,dvEdge,kiteAreasOnVertex,cellsOnVertex,nEdgesOnCell,areaCell)
end

function compute_weightsOnEdge_trisk!(weightsOnEdge,mesh::VoronoiMesh)
    nedges_max = max_length(eltype(mesh.edgesOnEdge))
    nedges = length(mesh.verticesOnEdge)
    size(weightsOnEdge) == (nedges_max,nedges) || throw(DimensionMismatch())
    return compute_weightsOnEdge_trisk!(weightsOnEdge,mesh.verticesOnEdge,mesh.cellsOnEdge,mesh.edgesOnEdge,mesh.dcEdge,mesh.dvEdge,mesh.kiteAreasOnVertex,mesh.cellsOnVertex,mesh.nEdgesOnCell,mesh.areaCell)
end

function compute_weightsOnEdge_trisk(mesh::VoronoiMesh)
    return compute_weightsOnEdge_trisk(mesh.verticesOnEdge,mesh.cellsOnEdge,mesh.edgesOnEdge,mesh.dcEdge,mesh.dvEdge,mesh.kiteAreasOnVertex,mesh.cellsOnVertex,mesh.nEdgesOnCell,mesh.areaCell)
end