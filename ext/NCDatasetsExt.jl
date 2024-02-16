module NCDatasetsExt

using VoronoiMeshDataStruct, NCDatasets, TensorsLite

function _on_a_sphere(ncfile::NCDatasets.NCDataset)
    oas = lowercase(ncfile.attrib["on_a_sphere"])
    if oas in ("yes","y")
        return Val(true)
    else
        return Val(false)
    end
end

function VoronoiMeshDataStruct.CellConnectivity(ncfile::NCDatasets.NCDataset)
    nEdgesOnCell = ncfile["nEdgesOnCell"][:]
    maxEdges = maximum(nEdgesOnCell)
    
    return CellConnectivity(maxEdges,ncfile)
end

function VoronoiMeshDataStruct.CellConnectivity(maxEdges::Integer,ncfile::NCDatasets.NCDataset)
    verticesOnCellArray = ncfile["verticesOnCell"]
    verticesOnCell = VariableLengthIndices{maxEdges}.((verticesOnCellArray[1:maxEdges,k]...,) for k in axes(verticesOnCellArray,2))

    edgesOnCellArray = ncfile["edgesOnCell"]
    edgesOnCell = VariableLengthIndices{maxEdges}.((edgesOnCellArray[1:maxEdges,k]...,) for k in axes(edgesOnCellArray,2))

    cellsOnCellArray = ncfile["cellsOnCell"]
    cellsOnCell = VariableLengthIndices{maxEdges}.((cellsOnCellArray[1:maxEdges,k]...,) for k in axes(cellsOnCellArray,2))

    return CellConnectivity(verticesOnCell,edgesOnCell,cellsOnCell)
end

function VoronoiMeshDataStruct.CellBase(ncfile::NCDatasets.NCDataset)
    nEdges = ncfile["nEdgesOnCell"][:]
    maxEdges = maximum(nEdges)
    indices = CellConnectivity(maxEdges,ncfile)
    
    position = VecArray(x=ncfile["xCell"][:],
                        y=ncfile["yCell"][:],
                        z=ncfile["zCell"][:])

    return CellBase(length(nEdges),indices,nEdges,position,_on_a_sphere(ncfile))
end

const cell_info_vectors = (longitude="lonCell", latitude="latCell",
                          meshDensity="meshDensity",indexToID="indexToCellID",
                          area="areaCell", bdyMask="bdyMaskCell")

const cell_info_matrices_max_edges = (defcA="defc_a", defcB="defc_b",
                                      xGradientCoeff="cell_gradient_coef_x", yGradientCoeff="cell_gradient_coef_y")

function VoronoiMeshDataStruct.CellInfo(ncfile::NCDatasets.NCDataset)
    cells = CellBase(ncfile)
    max_nedges = VoronoiMeshDataStruct.max_n_edges(typeof(cells))
    cellinfo = CellInfo(cells)

    for (field_name, nc_name) in pairs(cell_info_vectors)
        if haskey(ncfile,nc_name)
            setproperty!(cellinfo,field_name,ncfile[nc_name][:])
        end
    end

    if haskey(ncfile,"localVerticalUnitVectors")
        lvuva = ncfile["localVerticalUnitVectors"]
        cellinfo.verticalUnitVectors = VecArray(x=lvuva[1,:], y=lvuva[2,:], z=lvuva[3,:])
    end

    if haskey(ncfile,"cellTangentPlane")
        ctp = ncfile["cellTangentPlane"]
        cellinfo.tangentPlane = (VecArray(x=ctp[1,1,:],y=ctp[2,1,:],z=ctp[3,1,:]),VecArray(x=ctp[1,2,:],y=ctp[2,2,:],z=ctp[3,2,:]))
    end

    sl = Base.OneTo(max_nedges)
    for (field_name, nc_name) in pairs(cell_info_matrices_max_edges)
        if haskey(ncfile,nc_name)
            setproperty!(cellinfo,field_name,ncfile[nc_name][sl,:])
        end
    end

    if haskey(ncfile,"coeffs_reconstruct")
        coeffR = ncfile["coeffs_reconstruct"]
        sl = Base.OneTo(max_nedges)
        coeffsReconstruct = VecArray(x=coeffR[1,sl,:], y=coeffR[2,sl,:], z=coeffR[3,sl,:])
        cellinfo.coeffsReconstruct = coeffsReconstruct
    end

    return cellinfo    
end

function VoronoiMeshDataStruct.VertexConnectivity(ncfile::NCDatasets.NCDataset)
    edgesOnVertexArray = ncfile["edgesOnVertex"]
    edgesOnVertex = [(edgesOnVertexArray[1,k],edgesOnVertexArray[2,k],edgesOnVertexArray[3,k]) for k in axes(edgesOnVertexArray,2)]

    cellsOnVertexArray = ncfile["cellsOnVertex"]
    cellsOnVertex = [(cellsOnVertexArray[1,k],cellsOnVertexArray[2,k],cellsOnVertexArray[3,k]) for k in axes(cellsOnVertexArray,2)]
    return VertexConnectivity(edgesOnVertex,cellsOnVertex)
end

function VoronoiMeshDataStruct.VertexBase(ncfile::NCDatasets.NCDataset)
    indices = VertexConnectivity(ncfile)
    position = VecArray(x=ncfile["xVertex"][:],
                        y=ncfile["yVertex"][:],
                        z=ncfile["zVertex"][:])

    return VertexBase(length(position.x),indices,position,_on_a_sphere(ncfile))
end

const vertex_info_vectors = (longitude="lonVertex", latitude="latVertex",
                             indexToID="indexToVertexID",
                             area="areaTriangle", bdyMask="bdyMaskVertex")

function VoronoiMeshDataStruct.VertexInfo(ncfile::NCDatasets.NCDataset)
    vertexBase = VertexBase(ncfile)
    vertex = VertexInfo(vertexBase)

    for (field_name, nc_name) in pairs(vertex_info_vectors)
        if haskey(ncfile,nc_name)
            setproperty!(vertex,field_name,ncfile[nc_name][:])
        end
    end

    if haskey(ncfile,"kiteAreasOnVertex")
        kiteAreasArray = ncfile["kiteAreasOnVertex"]
        kiteAreas = [(kiteAreasArray[1,k],kiteAreasArray[2,k],kiteAreasArray[3,k]) for k in axes(kiteAreasArray,2)]
        vertex.kiteAreas = kiteAreas
    end

    return vertex
end

function VoronoiMeshDataStruct.EdgeConnectivity(ncfile::NCDatasets.NCDataset)
    verticesOnEdgeArray = ncfile["verticesOnEdge"]
    verticesOnEdge = [(verticesOnEdgeArray[1,k],verticesOnEdgeArray[2,k]) for k in axes(verticesOnEdgeArray,2)]

    cellsOnEdgeArray = ncfile["cellsOnEdge"]
    cellsOnEdge = [(cellsOnEdgeArray[1,k],cellsOnEdgeArray[2,k]) for k in axes(cellsOnEdgeArray,2)]
    return EdgeConnectivity(verticesOnEdge,cellsOnEdge)
end

function VoronoiMeshDataStruct.EdgeBase(ncfile::NCDatasets.NCDataset)
    indices = EdgeConnectivity(ncfile)
    position = VecArray(x=ncfile["xEdge"][:],
                        y=ncfile["yEdge"][:],
                        z=ncfile["zEdge"][:])

    return EdgeBase(length(position.x),indices,position,_on_a_sphere(ncfile))
end

function VoronoiMeshDataStruct.EdgeVelocityReconstruction(ncfile::NCDatasets.NCDataset)
    nEdges = ncfile["nEdgesOnEdge"][:]
    max_n_edges = maximum(nEdges)

    edgesOnEdgeArray = ncfile["edgesOnEdge"]
    indices = VariableLengthIndices{max_n_edges}.((edgesOnEdgeArray[1:max_n_edges,k]...,) for k in axes(edgesOnEdgeArray,2))
    weights = ncfile["weightsOnEdge"][Base.OneTo(max_n_edges),:]
    return EdgeVelocityReconstruction(nEdges,indices,weights)
end

const edge_info_vectors = (longitude="lonEdge", latitude="latEdge",
                           indexToID="indexToEdgeID", dv="dvEdge", dc="dcEdge",
                           angle="angleEdge", bdyMask="bdyMaskEdge")

function VoronoiMeshDataStruct.EdgeInfo(ncfile::NCDatasets.NCDataset)
    edgeinfo = EdgeInfo(EdgeBase(ncfile),EdgeVelocityReconstruction(ncfile))

    for (field_name, nc_name) in pairs(edge_info_vectors)
        if haskey(ncfile,nc_name)
            setproperty!(edgeinfo,field_name,ncfile[nc_name][:])
        end
    end

    if haskey(ncfile,"edgeNormalVectors")
        env = ncfile["edgeNormalVectors"]
        edgeinfo.normalVectors = VecArray(x=env[1,:], y=env[2,:], z=env[3,:])
    end

    if haskey(ncfile,"deriv_two")
        edgeinfo.derivTwo = ncfile["deriv_two"][:,:,:]
    end

    return edgeinfo
end

function VoronoiMeshDataStruct.VoronoiMesh(ncfile::NCDatasets.NCDataset)
    return VoronoiMesh(VertexInfo(ncfile),CellInfo(ncfile),EdgeInfo(ncfile)) 
end

end