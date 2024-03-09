module VoronoiMeshDataStruct

using TensorsLite

export VariableLengthIndices
export CellConnectivity, CellBase, CellInfo
export VertexConnectivity, VertexBase, VertexInfo
export EdgeConnectivity, EdgeBase, EdgeVelocityReconstruction, EdgeInfo
export VoronoiMesh

export on_a_sphere, max_n_edges, max_n_edges_vel_reconstruction, float_precision, integer_precision
export create_cells_polygons, create_cells_polygons_periodic

include("variable_length_indices.jl")

include("cell_datastruct.jl")

include("vertices_datastruct.jl")

include("edges_datastruct.jl")

struct VoronoiMesh{S,N,N2,TI,TF,Tz}
    vertices::VertexInfo{S,TI,TF,Tz}
    cells::CellInfo{S,N,TI,TF,Tz}
    edges::EdgeInfo{S,N2,TI,TF,Tz}
    attributes::Dict{Symbol,Union{String,Float64}}
end

Base.getproperty(mesh::VoronoiMesh,s::Symbol) = _getproperty(mesh,Val(s))
_getproperty(mesh::VoronoiMesh,::Val{s}) where s = getfield(mesh,s)

for s in (:nVertices,:cellsOnVertex,:edgesOnVertex,:xVertex,:yVertex,:zVertex, :lonVertex, :latVertex,
          :indexToVertexID, :areaTriangle, :bdyMaskVertex, :kiteAreasOnVertex)
    @eval _getproperty(mesh::VoronoiMesh,::Val{$(QuoteNode(s))}) = getproperty(getfield(mesh,:vertices),$(QuoteNode(s)))
end

for s in (:nCells,:verticesOnCell,:edgesOnCell,:cellsOnCell,:xCell,:yCell,:zCell,:nEdgesOnCell,
          :lonCell, :latCell, :meshDensity, :indexToCellID, :areaCell, :bdyMaskCell,
          :localVerticalUnitVectors, :cellTangentPlane, :defc_a,
          :defc_b, :cell_gradient_coef_x, :cell_gradient_coef_y, :coeffs_reconstruct)
    @eval _getproperty(mesh::VoronoiMesh,::Val{$(QuoteNode(s))}) = getproperty(getfield(mesh,:cells),$(QuoteNode(s)))
end

for s in (:nEdges,:verticesOnEdge,:cellsOnEdge,:xEdge,:yEdge,:zEdge, :nEdgesOnEdge, :edgesOnEdge, :weightsOnEdge,:lonEdge,
          :latEdge, :indexToEdgeID, :dvEdge, :bdyMaskEdge, :dcEdge, :angleEdge, :edgeNormalVectors, :deriv_two)
    @eval _getproperty(mesh::VoronoiMesh,::Val{$(QuoteNode(s))}) = getproperty(getfield(mesh,:edges),$(QuoteNode(s)))
end

function create_cells_polygons end
function create_cells_polygons_periodic end

end
