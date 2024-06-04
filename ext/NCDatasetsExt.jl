module NCDatasetsExt

using VoronoiMeshDataStruct, NCDatasets, TensorsLite, Zeros, ImmutableVectors

function _on_a_sphere(ncfile::NCDatasets.NCDataset)
    oas = lowercase(strip(ncfile.attrib["on_a_sphere"]))
    if oas in ("yes","y")
        return (Val(true), true)
    else
        return (Val(false), false)
    end
end

function construct_elementsOnElement(n::Val{maxEdges},elementsOnElementArray::AbstractMatrix{TE},nElemtensOnElement::AbstractVector{TI}) where {maxEdges,TE,TI}
    elementsOnElement = Vector{ImmutableVector{maxEdges,TE}}(undef,size(elementsOnElementArray,2))
    @inbounds for k in axes(elementsOnElementArray,2)
        elementsOnElement[k] = ImmutableVector{maxEdges}(ntuple(i->(@inbounds elementsOnElementArray[i,k]),n),nElemtensOnElement[k])
    end
    return elementsOnElement
end

function VoronoiMeshDataStruct.CellConnectivity(me::Val{maxEdges},nEdgesOnCell::AbstractVector{TI},ncfile::NCDatasets.NCDataset) where {maxEdges,TI}
    verticesOnCellArray = ncfile["verticesOnCell"][:,:]
    verticesOnCell = construct_elementsOnElement(me,verticesOnCellArray,nEdgesOnCell)

    edgesOnCellArray = ncfile["edgesOnCell"][:,:]
    edgesOnCell = construct_elementsOnElement(me,edgesOnCellArray,nEdgesOnCell)

    cellsOnCellArray = ncfile["cellsOnCell"][:,:]
    cellsOnCell = construct_elementsOnElement(me,cellsOnCellArray,nEdgesOnCell)

    return CellConnectivity(verticesOnCell,edgesOnCell,cellsOnCell)
end

for N in 6:12
    for T in (Int64,Int32)
        for TE in (Int64,Int32,Float64,Float32)
            precompile(construct_elementsOnElement,(Val{N},Matrix{TE},Vector{T}))
        end
        precompile(VoronoiMeshDataStruct.CellConnectivity,(Val{N},Vector{T},NCDatasets.NCDataset{Nothing,Missing}))
    end
end

function VoronoiMeshDataStruct.CellConnectivity(ncfile::NCDatasets.NCDataset)
    nEdgesOnCell = ncfile["nEdgesOnCell"][:]::Vector{Int32}
    maxEdges = Int(maximum(nEdgesOnCell))
    
    return CellConnectivity(Val(maxEdges),nEdgesOnCell,ncfile)
end

precompile(VoronoiMeshDataStruct.CellConnectivity,(NCDatasets.NCDataset{Nothing,Missing},))

function VoronoiMeshDataStruct.CellBase(onSphere::Val{on_sphere},mE::Val{maxEdges},nEdges,ncfile::NCDatasets.NCDataset) where {on_sphere,maxEdges}
    indices = CellConnectivity(mE,nEdges,ncfile)
    position = on_sphere ? VecArray(x=ncfile["xCell"][:], y=ncfile["yCell"][:], z=ncfile["zCell"][:]) : VecArray(x=ncfile["xCell"][:], y=ncfile["yCell"][:])

    return CellBase(length(nEdges),indices,nEdges,position,onSphere)
end

for N in 6:12
    for on_sphere in (true,false)
        precompile(VoronoiMeshDataStruct.CellBase,(Val{on_sphere}, Val{N}, Vector{Int32}, NCDatasets.NCDataset{Nothing,Missing}))
    end
end

function VoronoiMeshDataStruct.CellBase(ncfile::NCDatasets.NCDataset)
    nEdges = ncfile["nEdgesOnCell"][:]::Vector{Int32}
    maxEdges = Int(maximum(nEdges))
    onSphere, _ = _on_a_sphere(ncfile)
    return CellBase(onSphere,Val(maxEdges),nEdges,ncfile)
end
precompile(VoronoiMeshDataStruct.CellBase,(NCDatasets.NCDataset{Nothing,Missing},))

const cell_info_vectors = (longitude="lonCell", latitude="latCell",
                          meshDensity="meshDensity",indexToID="indexToCellID",
                          area="areaCell", bdyMask="bdyMaskCell")

const cell_info_matrices_max_edges = (defcA="defc_a", defcB="defc_b",
                                      xGradientCoeff="cell_gradient_coef_x", yGradientCoeff="cell_gradient_coef_y")

function VoronoiMeshDataStruct.CellInfo(ncfile::NCDatasets.NCDataset)
    cells = CellBase(ncfile)
    return CellInfo(cells,ncfile)
end
precompile(VoronoiMeshDataStruct.CellInfo,(NCDatasets.NCDataset{Nothing,Missing},))

function VoronoiMeshDataStruct.CellInfo(cells::CellBase{on_sphere,max_nedges},ncfile::NCDatasets.NCDataset) where {on_sphere,max_nedges}

    cellinfo = CellInfo(cells)

    for (field_name, nc_name) in pairs(cell_info_vectors)
        if haskey(ncfile,nc_name)
            setproperty!(cellinfo,field_name,ncfile[nc_name][:])
        end
    end

    if haskey(ncfile,"localVerticalUnitVectors")
        lvuva = ncfile["localVerticalUnitVectors"]
        cellinfo.verticalUnitVectors = on_sphere ? VecArray(x=lvuva[1,:], y=lvuva[2,:], z=lvuva[3,:]) : VecArray(z=lvuva[3,:])
    end

    if haskey(ncfile,"cellTangentPlane")
        ctp = ncfile["cellTangentPlane"]
        cellinfo.tangentPlane = on_sphere ?
                                (VecArray(x=ctp[1,1,:],y=ctp[2,1,:],z=ctp[3,1,:]),VecArray(x=ctp[1,2,:],y=ctp[2,2,:],z=ctp[3,2,:])) :
                                (VecArray(x=ctp[1,1,:],y=ctp[2,1,:]),VecArray(x=ctp[1,2,:],y=ctp[2,2,:]))
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
        coeffsReconstruct = on_sphere ?
                            VecArray(x=coeffR[1,sl,:], y=coeffR[2,sl,:], z=coeffR[3,sl,:]) :
                            VecArray(x=coeffR[1,sl,:], y=coeffR[2,sl,:])
        cellinfo.coeffsReconstruct = coeffsReconstruct
    end

    return cellinfo    
end

for N in 6:12
    for on_sphere in (true,false)
        for TF in (Float64,Float32)
            for Tz in (Float64,Float32,Zero)
                precompile(VoronoiMeshDataStruct.CellInfo,(VoronoiMeshDataStruct.CellBase{on_sphere,N,Int32,TF,Tz}, NCDatasets.NCDataset{Nothing,Missing}))
            end
        end
    end
end

function VoronoiMeshDataStruct.VertexConnectivity(ncfile::NCDatasets.NCDataset)
    edgesOnVertexArray = ncfile["edgesOnVertex"][:,:]::Matrix{Int32}
    edgesOnVertex = [(edgesOnVertexArray[1,k],edgesOnVertexArray[2,k],edgesOnVertexArray[3,k]) for k in axes(edgesOnVertexArray,2)]

    cellsOnVertexArray = ncfile["cellsOnVertex"][:,:]::Matrix{Int32}
    cellsOnVertex = [(cellsOnVertexArray[1,k],cellsOnVertexArray[2,k],cellsOnVertexArray[3,k]) for k in axes(cellsOnVertexArray,2)]
    return VertexConnectivity(edgesOnVertex,cellsOnVertex)
end

precompile(VoronoiMeshDataStruct.VertexConnectivity,(NCDatasets.NCDataset{Nothing,Missing},))

function VoronoiMeshDataStruct.VertexBase(ncfile::NCDatasets.NCDataset)
    indices = VertexConnectivity(ncfile)

    onSphere, on_sphere = _on_a_sphere(ncfile)

    position = on_sphere ? VecArray(x=ncfile["xVertex"][:],y=ncfile["yVertex"][:],z=ncfile["zVertex"][:]) : VecArray(x=ncfile["xVertex"][:],y=ncfile["yVertex"][:])

    return VertexBase(length(position.x),indices,position,onSphere)
end

precompile(VoronoiMeshDataStruct.VertexBase,(NCDatasets.NCDataset{Nothing,Missing},))

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

precompile(VoronoiMeshDataStruct.VertexInfo,(NCDatasets.NCDataset{Nothing,Missing},))

function VoronoiMeshDataStruct.EdgeConnectivity(ncfile::NCDatasets.NCDataset)
    verticesOnEdgeArray = ncfile["verticesOnEdge"][:,:]::Matrix{Int32}
    verticesOnEdge = [(verticesOnEdgeArray[1,k],verticesOnEdgeArray[2,k]) for k in axes(verticesOnEdgeArray,2)]

    cellsOnEdgeArray = ncfile["cellsOnEdge"][:,:]::Matrix{Int32}
    cellsOnEdge = [(cellsOnEdgeArray[1,k],cellsOnEdgeArray[2,k]) for k in axes(cellsOnEdgeArray,2)]
    return EdgeConnectivity(verticesOnEdge,cellsOnEdge)
end

precompile(VoronoiMeshDataStruct.EdgeConnectivity,(NCDatasets.NCDataset{Nothing,Missing},))

function VoronoiMeshDataStruct.EdgeBase(ncfile::NCDatasets.NCDataset)
    indices = EdgeConnectivity(ncfile)

    onSphere, on_sphere = _on_a_sphere(ncfile)

    position = on_sphere ? VecArray(x=ncfile["xEdge"][:],y=ncfile["yEdge"][:],z=ncfile["zEdge"][:]) : VecArray(x=ncfile["xEdge"][:],y=ncfile["yEdge"][:])

    return EdgeBase(length(position.x),indices,position,onSphere)
end

precompile(VoronoiMeshDataStruct.EdgeBase,(NCDatasets.NCDataset{Nothing,Missing},))

function VoronoiMeshDataStruct.EdgeVelocityReconstruction(ncfile::NCDatasets.NCDataset)
    nEdges = ncfile["nEdgesOnEdge"][:]::Vector{Int32}
    max_n_edges = Int(maximum(nEdges))
    return EdgeVelocityReconstruction(Val(max_n_edges),nEdges,ncfile)
end

precompile(VoronoiMeshDataStruct.EdgeVelocityReconstruction,(NCDatasets.NCDataset{Nothing,Missing},))

function VoronoiMeshDataStruct.EdgeVelocityReconstruction(ne::Val{max_n_edges},nEdges,ncfile::NCDatasets.NCDataset) where {max_n_edges}
    edgesOnEdgeArray = ncfile["edgesOnEdge"][:,:]::Matrix{Int32}
    indices = construct_elementsOnElement(ne,edgesOnEdgeArray,nEdges)
    weightsArray = ncfile["weightsOnEdge"][:,:]
    weights = construct_elementsOnElement(ne,weightsArray,nEdges)
    return EdgeVelocityReconstruction(nEdges,indices,weights)
end

for N in 10:24
    precompile(VoronoiMeshDataStruct.EdgeVelocityReconstruction,(Val{N},Vector{Int32},NCDatasets.NCDataset{Nothing,Missing}))
end

const edge_info_vectors = (longitude="lonEdge", latitude="latEdge",
                           indexToID="indexToEdgeID", dv="dvEdge", dc="dcEdge",
                           angle="angleEdge", bdyMask="bdyMaskEdge")

function VoronoiMeshDataStruct.EdgeInfo(ncfile::NCDatasets.NCDataset)
    edgeinfo = EdgeInfo(EdgeBase(ncfile),EdgeVelocityReconstruction(ncfile))

    _, on_sphere = _on_a_sphere(ncfile)

    for (field_name, nc_name) in pairs(edge_info_vectors)
        if haskey(ncfile,nc_name)
            setproperty!(edgeinfo,field_name,ncfile[nc_name][:])
        end
    end

    if haskey(ncfile,"edgeNormalVectors")
        env = ncfile["edgeNormalVectors"]
        edgeinfo.normalVectors = on_sphere ? VecArray(x=env[1,:], y=env[2,:], z=env[3,:]) : VecArray(x=env[1,:], y=env[2,:])
    end

    if haskey(ncfile,"deriv_two")
        edgeinfo.derivTwo = ncfile["deriv_two"][:,:,:]
    end

    return edgeinfo
end

precompile(VoronoiMeshDataStruct.EdgeInfo,(NCDatasets.NCDataset{Nothing,Missing},))

function VoronoiMeshDataStruct.VoronoiMesh(ncfile::NCDatasets.NCDataset)
    attributes = Dict{Symbol,Union{String,Float64,Float32,Int64,Int32}}()
    for (key,val) in ncfile.attrib
        val isa String && (val = String(strip(val)))
        attributes[Symbol(key)] = val
    end
    return VoronoiMesh(VertexInfo(ncfile),CellInfo(ncfile),EdgeInfo(ncfile),attributes) 
end

precompile(VoronoiMeshDataStruct.VoronoiMesh,(NCDatasets.NCDataset{Nothing,Missing},))

for func in  (:VoronoiMesh,
              :CellInfo,:CellBase,:CellConnectivity,
              :VertexInfo,:VertexBase,:VertexConnectivity,
              :EdgeInfo,:EdgeBase,:EdgeConnectivity,:EdgeVelocityReconstruction)

    @eval begin
        VoronoiMeshDataStruct.$func(file_name::String) = NCDataset(file_name) do f; VoronoiMeshDataStruct.$func(f);end
        precompile(VoronoiMeshDataStruct.$func,(String,))
    end
end

end