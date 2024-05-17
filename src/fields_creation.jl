@inline function check_sizes(size_expected,this_size)
    this_size == size_expected || throw(DimensionMismatch("Output array has incorrect size.\nExpected size: $size_expected\nGot: $this_size"))
end

function compute_edge_position_periodic!(epos,cpos,cellsOnEdge,xp::Number,yp::Number)
    check_sizes(size(cellsOnEdge),size(epos))

    @inbounds for i in eachindex(cellsOnEdge)
        ic1,ic2 = cellsOnEdge[i]
        cpos2 = cpos[ic2]
        cpos1 = closest(cpos2,cpos[ic1],xp,yp)
        epos[i] = 0.5*(cpos2 + cpos1)
    end
    return epos
end

compute_edge_position!(epos,edges::EdgeBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_edge_position_periodic!(epos,cells.position,edges.indices.cells,xp,yp)

compute_edge_position!(epos,mesh::VoronoiMesh{false}) = compute_edge_position!(epos,mesh.edges.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_edge_position_periodic(cpos,cellsOnEdge,xp::Number,yp::Number)
    epos = similar(cpos,size(cellsOnEdge))
    return compute_edge_position_periodic!(epos,cpos,cellsOnEdge,xp,yp)
end

compute_edge_position(edges::EdgeBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_edge_position_periodic(cells.position,edges.indices.cells,xp,yp)
compute_edge_position!(edges::EdgeBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_edge_position_periodic!(edges.position,cells.position,edges.indices.cells,xp,yp)

compute_edge_position(mesh::VoronoiMesh{false}) = compute_edge_position(mesh.edges.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)
compute_edge_position!(mesh::VoronoiMesh) = compute_edge_position!(mesh.edges.position,mesh)

function compute_edge_normals_periodic!(n,cpos,cellsOnEdge,xp::Number,yp::Number)
    check_sizes(size(cellsOnEdge),size(n))

    @inbounds for i in eachindex(cellsOnEdge)
        ic1,ic2 = cellsOnEdge[i]
        cpos2 = cpos[ic2]
        cpos1 = closest(cpos2,cpos[ic1],xp,yp)
        n[i] = normalize(cpos2 - cpos1)
    end

    return n
end

compute_edge_normals!(n,edges::EdgeBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_edge_normals_periodic!(n,cells.position,edges.indices.cells,xp,yp)

compute_edge_normals!(n,mesh::VoronoiMesh{false}) = compute_edge_normals!(n,mesh.edges.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

#Since the edges are always equidistant to both cell positions, the vector connecting them will be parallel to the arc tangent at the edge, right??
function compute_edge_normals_on_sphere!(n,cpos,cellsOnEdge)
    check_sizes(size(cellsOnEdge),size(n))

    @inbounds for i in eachindex(cellsOnEdge)
        ic1,ic2 = cellsOnEdge[i]
        n[i] = normalize(cpos[ic2] - cpos[ic1])
    end
    return n
end

compute_edge_normals!(n,edges::EdgeBase{true},cells::CellBase{true}) = compute_edge_normals_on_sphere!(n,cells.position,edges.indices.cells)

compute_edge_normals!(n,mesh::VoronoiMesh{true}) = compute_edge_normals!(n,mesh.edges.base,mesh.cells.base)

function compute_edge_normals_periodic(cpos,cellsOnEdge,xp::Number,yp::Number)
    n = similar(cpos,size(cellsOnEdge))
    return compute_edge_normals_periodic!(n,cpos,cellsOnEdge,xp,yp)
end

compute_edge_normals(edges::EdgeBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_edge_normals_periodic(cells.position,edges.indices.cells,xp,yp)

compute_edge_normals(mesh::VoronoiMesh{false}) = compute_edge_normals(mesh.edges.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_edge_normals_on_sphere(cpos,cellsOnEdge)
    n = similar(cpos,size(cellsOnEdge))
    return compute_edge_normals_on_sphere!(n,cpos,cellsOnEdge)
end

compute_edge_normals(edges::EdgeBase{true},cells::CellBase{true}) = compute_edge_normals_on_sphere(cells.position,edges.indices.cells)
compute_edge_normals(mesh::VoronoiMesh{true}) = compute_edge_normals_on_sphere(mesh.cells.position,mesh.edges.indices.cells)

function compute_edge_normals!(mesh::VoronoiMesh)
    if isdefined(mesh.edge,:normalVectors)
        compute_edge_normals!(mesh.edge.normalVectors,mesh)
    else
        mesh.edge.normalVectors = compute_edge_normals(mesh)
    end
    return mesh.edge.normalVectors
end

function compute_area_triangles_periodic!(output,cpos,cellsOnVertex,xp::Number,yp::Number)
    check_sizes(size(cellsOnVertex),size(output))

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

compute_area_triangles(vertices::VertexBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_area_triangles_periodic(cells.position,vertices.indices.cells,xp,yp)

compute_area_triangles(mesh::VoronoiMesh{false}) = compute_area_triangles(mesh.vertices.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

compute_area_triangles!(output,vertices::VertexBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_area_triangles_periodic!(output,cells.position,vertices.indices.cells,xp,yp)

compute_area_triangles!(output,mesh::VoronoiMesh{false}) = compute_area_triangles!(output,mesh.vertices.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_area_triangles!(mesh::VoronoiMesh)
    if isdefined(mesh.vertices,:area)
        compute_area_triangles!(mesh.vertices.area,mesh)
    else
        mesh.vertices.area = compute_area_triangles(mesh)
    end
    return mesh.vertices.area
end

function compute_kite_areas_periodic!(output,cpos,vpos,cellsOnVertex,xp,yp)
    check_sizes(size(cellsOnVertex),size(output))

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

compute_kite_areas!(output,vertices::VertexBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_kite_areas_periodic!(output,cells.position,vertices.position,vertices.indices.cells,xp,yp)

compute_kite_areas!(output,mesh::VoronoiMesh{false}) = compute_kite_areas_periodic!(output,mesh.cells.position,mesh.vertices.position,mesh.vertices.indices.cells,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_kite_areas_periodic(cpos,vpos,cellsOnVertex,xp,yp)
    output = Vector{NTuple{3,nonzero_eltype(eltype(cpos))}}(undef,length(cellsOnVertex))
    return compute_kite_areas_periodic!(output,cpos,vpos,cellsOnVertex,xp,yp)
end

compute_kite_areas(vertices::VertexBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_kite_areas_periodic(cells.position,vertices.position,vertices.indices.cells,xp,yp)

compute_kite_areas(mesh::VoronoiMesh{false}) = compute_kite_areas(mesh.vertices.base, mesh.cells.base, mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_area_cells_periodic!(output,vpos,verticesOnCell,xp::Number,yp::Number)
    check_sizes(size(verticesOnCell),size(output))

    @inbounds for c in eachindex(verticesOnCell)
        output[c] = area(vpos,verticesOnCell[c],xp,yp)
    end
    return output
end

compute_area_cells!(output,cells::CellBase{false},vertices::VertexBase{false},xp::Number,yp::Number) = compute_area_cells_periodic!(output,vertices.position,cells.indices.vertices,xp,yp)

compute_area_cells!(output,mesh::VoronoiMesh{false}) = compute_area_cells!(output,mesh.cells.base,mesh.vertices.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_area_cells_periodic(vpos,verticesOnCell,xp::Number,yp::Number)
    output = Vector{nonzero_eltype(eltype(vpos))}(undef,length(verticesOnCell))
    return compute_area_cells_periodic!(output,vpos,verticesOnCell,xp,yp)
end

compute_area_cells(cells::CellBase{false},vertices::VertexBase{false},xp::Number,yp::Number) = compute_area_cells_periodic(vertices.position,cells.indices.vertices,xp,yp)

compute_area_cells(mesh::VoronoiMesh{false}) = compute_area_cells(mesh.cells.base,mesh.vertices.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_area_cells!(mesh::VoronoiMesh)
    if isdefined(mesh.cells,:area)
        compute_area_cells!(mesh.cells.area,mesh)
    else
        mesh.cells.area = compute_area_cells(mesh)
    end
    return mesh.cells.area
end

function compute_dcEdge_periodic!(output,cpos,cellsOnEdge,xp::Number,yp::Number)
    check_sizes(size(cellsOnEdge),size(output))

    @inbounds for e in eachindex(cellsOnEdge)
        c1,c2 = cellsOnEdge[e]
        c1pos = cpos[c1]
        c2pos = closest(c1pos,cpos[c2],xp,yp)
        output[e] = norm(c2pos-c1pos)
    end
    return output
end

compute_dcEdge!(output,edges::EdgeBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_dcEdge_periodic!(output,cells.position,edges.indices.cells,xp,yp)

compute_dcEdge!(output,mesh::VoronoiMesh{false}) = compute_dcEdge!(output,mesh.edges.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_dcEdge_periodic(cpos,cellsOnEdge,xp::Number,yp::Number)
    output = Vector{nonzero_eltype(eltype(cpos))}(undef,length(cellsOnEdge))
    return compute_dcEdge_periodic!(output,cpos,cellsOnEdge,xp,yp)
end

compute_dcEdge(edges::EdgeBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_dcEdge_periodic(cells.position,edges.indices.cells,xp,yp)

compute_dcEdge(mesh::VoronoiMesh{false}) = compute_dcEdge(mesh.edges.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_dcEdge!(mesh::VoronoiMesh) 
    if isdefined(mesh.edges,:dc)
        return compute_dcEdge!(mesh.edges.dc,mesh)
    else
        mesh.edges.dc = compute_dcEdge(mesh)
        return mesh.edges.dc
    end
end

function compute_dvEdge_periodic!(output,vpos,verticesOnEdge,xp::Number,yp::Number)
    check_sizes(size(verticesOnEdge),size(output))

    @inbounds for e in eachindex(verticesOnEdge)
        v1,v2 = verticesOnEdge[e]
        v1pos = vpos[v1]
        v2pos = closest(v1pos,vpos[v2],xp,yp)
        output[e] = norm(v2pos-v1pos)
    end
    return output
end

compute_dvEdge!(output,edges::EdgeBase{false},vertices::VertexBase{false},xp::Number,yp::Number) = compute_dvEdge_periodic!(output,vertices.position,edges.indices.vertices,xp,yp)

compute_dvEdge!(output,mesh::VoronoiMesh{false}) = compute_dvEdge!(output,mesh.edges.base,mesh.vertices.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_dvEdge_periodic(vpos,verticesOnEdge,xp::Number,yp::Number)
    output = Vector{nonzero_eltype(eltype(vpos))}(undef,length(verticesOnEdge))
    return compute_dvEdge_periodic!(output,vpos,verticesOnEdge,xp,yp)
end

compute_dvEdge(edges::EdgeBase{false},vertices::VertexBase{false},xp::Number,yp::Number) = compute_dvEdge_periodic(vertices.position,edges.indices.vertices,xp,yp)

compute_dvEdge(mesh::VoronoiMesh{false}) = compute_dvEdge(mesh.edges.base,mesh.vertices.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_dvEdge!(mesh::VoronoiMesh) 
    if isdefined(mesh.edges,:dv)
        return compute_dvEdge!(mesh.edges.dv,mesh)
    else
        mesh.edges.dv = compute_dvEdge(mesh)
        return mesh.edges.dv
    end
end

function compute_angleEdge_periodic!(output,cpos,cellsOnEdge,xp::Number,yp::Number)
    check_sizes(size(cellsOnEdge),size(output))

    @inbounds for i in eachindex(cellsOnEdge)
        ic1,ic2 = cellsOnEdge[i]
        cpos2 = cpos[ic2]
        cpos1 = closest(cpos2,cpos[ic1],xp,yp)
        output[i] = acos(normalize(cpos2 - cpos1) ⋅ 𝐢)
    end
    return output
end

compute_angleEdge!(output,edges::EdgeBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_angleEdge_periodic!(output,cells.position,edges.indices.cells,xp,yp)

compute_angleEdge!(output,mesh::VoronoiMesh{false}) = compute_angleEdge!(output,mesh.edges.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

function compute_angleEdge_periodic(cpos,cellsOnEdge,xp::Number,yp::Number)
    output = Vector{nonzero_eltype(eltype(cpos))}(undef,length(cellsOnEdge))
    return compute_angleEdge_periodic!(output,cpos,cellsOnEdge,xp,yp)
end

compute_angleEdge(edges::EdgeBase{false},cells::CellBase{false},xp::Number,yp::Number) = compute_angleEdge_periodic(cells.position,edges.indices.cells,xp,yp)

compute_angleEdge(mesh::VoronoiMesh{false}) = compute_angleEdge(mesh.edges.base,mesh.cells.base,mesh.attributes[:x_period]::Float64,mesh.attributes[:y_period]::Float64)

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

"""
    select_kite_area(kiteAreaOnVertex::Vector{NTuple{3}}, cellsOnVertex::Vector{NTuple{3,<:Integer}}, v_i::Integer, c_i::Integer)
Returns the kite Area associated with cell `c_i` and vertex `v_i`.
Throws an error if vertex `v_i` doens't belong to cell `c_i`.
"""
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
    check_sizes((max_length(eltype(edgesOnEdge)),length(edgesOnEdge)),size(weightsOnEdge))

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

function compute_weightsOnEdge_trisk!(w,edges::EdgeBase,velRecon::EdgeVelocityReconstruction,cells::CellBase,vertices::VertexBase,dcEdge,dvEdge,kiteAreasOnVertex,areaCell)
    verticesOnEdge = edges.indices.vertices
    cellsOnEdge = edges.indices.cells
    edgesOnEdge = velRecon.indices
    cellsOnVertex = vertices.indices.cells
    nEdgesOnCell = cells.nEdges

    return compute_weightsOnEdge_trisk!(w,verticesOnEdge,cellsOnEdge,edgesOnEdge,dcEdge,dvEdge,kiteAreasOnVertex,cellsOnVertex,nEdgesOnCell,areaCell)
end

compute_weightsOnEdge_trisk!(weightsOnEdge,mesh::VoronoiMesh) = compute_weightsOnEdge_trisk!(weightsOnEdge,mesh.edges.base,mesh.edges.velRecon,mesh.cells.base,mesh.vertices.base,mesh.dcEdge,mesh.dvEdge,mesh.kiteAreasOnVertex,mesh.areaCell)

function compute_weightsOnEdge_trisk(verticesOnEdge,cellsOnEdge,edgesOnEdge,dcEdge,dvEdge,kiteAreasOnVertex,cellsOnVertex,nEdgesOnCell,areaCell)
    T = eltype(eltype(kiteAreasOnVertex))
    nedges_max = max_length(eltype(edgesOnEdge))
    nedges = length(verticesOnEdge)
    weightsOnEdge = zeros(T,nedges_max,nedges)
    return compute_weightsOnEdge_trisk!(weightsOnEdge,verticesOnEdge,cellsOnEdge,edgesOnEdge,dcEdge,dvEdge,kiteAreasOnVertex,cellsOnVertex,nEdgesOnCell,areaCell)
end

function compute_weightsOnEdge_trisk(mesh::VoronoiMesh)
    return compute_weightsOnEdge_trisk(mesh.verticesOnEdge,mesh.cellsOnEdge,mesh.edgesOnEdge,mesh.dcEdge,mesh.dvEdge,mesh.kiteAreasOnVertex,mesh.cellsOnVertex,mesh.nEdgesOnCell,mesh.areaCell)
end
