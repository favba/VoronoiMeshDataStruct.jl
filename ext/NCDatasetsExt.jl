module NCDatasetsExt

using VoronoiMeshDataStruct, NCDatasets, TensorsLite

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

function VoronoiMeshDataStruct.Cells(ncfile::NCDatasets.NCDataset)
    nEdges = ncfile["nEdgesOnCell"][:]
    maxEdges = maximum(nEdges)
    indices = CellConnectivity(maxEdges,ncfile)
    
    position = VecArray(x=ncfile["xCell"][:],
                        y=ncfile["yCell"][:],
                        z=ncfile["zCell"][:])

    return Cells(indices,nEdges,position)
end

function VoronoiMeshDataStruct.VertexConnectivity(ncfile::NCDatasets.NCDataset)
    edgesOnVertexArray = ncfile["edgesOnVertex"]
    edgesOnVertex = [(edgesOnVertexArray[1,k],edgesOnVertexArray[2,k],edgesOnVertexArray[3,k]) for k in axes(edgesOnVertexArray,2)]

    cellsOnVertexArray = ncfile["cellsOnVertex"]
    cellsOnVertex = [(cellsOnVertexArray[1,k],cellsOnVertexArray[2,k],cellsOnVertexArray[3,k]) for k in axes(cellsOnVertexArray,2)]
    return VertexConnectivity(edgesOnVertex,cellsOnVertex)
end

function VoronoiMeshDataStruct.Vertices(ncfile::NCDatasets.NCDataset)
    indices = VertexConnectivity(ncfile)
    position = VecArray(x=ncfile["xVertex"][:],
                        y=ncfile["yVertex"][:],
                        z=ncfile["zVertex"][:])

    return Vertices(indices,position)
end

function VoronoiMeshDataStruct.EdgeConnectivity(ncfile::NCDatasets.NCDataset)
    verticesOnEdgeArray = ncfile["verticesOnEdge"]
    verticesOnEdge = [(verticesOnEdgeArray[1,k],verticesOnEdgeArray[2,k]) for k in axes(verticesOnEdgeArray,2)]

    cellsOnEdgeArray = ncfile["cellsOnEdge"]
    cellsOnEdge = [(cellsOnEdgeArray[1,k],cellsOnEdgeArray[2,k]) for k in axes(cellsOnEdgeArray,2)]
    return EdgeConnectivity(verticesOnEdge,cellsOnEdge)
end

function VoronoiMeshDataStruct.Edges(ncfile::NCDatasets.NCDataset)
    indices = EdgeConnectivity(ncfile)
    position = VecArray(x=ncfile["xEdge"][:],
                        y=ncfile["yEdge"][:],
                        z=ncfile["zEdge"][:])

    return Edges(indices,position)
end


end