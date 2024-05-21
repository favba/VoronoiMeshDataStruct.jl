module VoronoiMeshDataStruct

using TensorsLite, TensorsLiteGeometry

export VariableLengthIndices
export CellConnectivity, CellBase, CellInfo
export VertexConnectivity, VertexBase, VertexInfo
export EdgeConnectivity, EdgeBase, EdgeVelocityReconstruction, EdgeInfo
export VoronoiMesh

export on_a_sphere, max_n_edges, max_n_edges_vel_reconstruction, float_precision, integer_precision

export compute_edge_position, compute_edge_position!, compute_edge_position_periodic, compute_edge_position_periodic!
export compute_edge_normals, compute_edge_normals!, compute_edge_normals_periodic, compute_edge_normals_periodic!
export compute_area_triangles, compute_area_triangles!, compute_kite_areas, compute_kite_areas!, compute_area_cells, compute_area_cells!
export compute_area_triangles_periodic, compute_area_triangles_periodic!
export compute_area_cells_periodic, compute_area_cells_periodic!
export compute_kite_areas_periodic, compute_kite_areas_periodic!
export compute_dcEdge, compute_dcEdge!
export compute_dcEdge_periodic, compute_dcEdge_periodic!
export compute_dvEdge, compute_dvEdge!
export compute_dvEdge_periodic, compute_dvEdge_periodic!
export compute_angleEdge, compute_angleEdge!
export compute_angleEdge_periodic, compute_angleEdge_periodic!
export compute_weightsOnEdge_trisk, compute_weightsOnEdge_trisk!

export create_cells_polygons, create_cells_polygons_periodic, create_dual_triangles, create_dual_triangles_periodic, create_edge_quadrilaterals, create_edge_quadrilaterals_periodic
export create_cell_linesegments, create_cell_linesegments_periodic

export graph_partition, find_obtuse_triangles
export select_kite_area, periodic_edges_mask

include("variable_length_indices.jl")

include("cell_datastruct.jl")

include("vertices_datastruct.jl")

include("edges_datastruct.jl")

struct VoronoiMesh{S,N,N2,TI,TF,Tz}
    vertices::VertexInfo{S,TI,TF,Tz}
    cells::CellInfo{S,N,TI,TF,Tz}
    edges::EdgeInfo{S,N2,TI,TF,Tz}
    attributes::Dict{Symbol,Union{String,Float64,Float32,Int64,Int32}}
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

include("fields_creation.jl")
include("utils.jl")

function create_cells_polygons end
function create_cells_polygons_periodic end

function create_dual_triangles end
function create_dual_triangles_periodic end

function create_edge_quadrilaterals end
function create_edge_quadrilaterals_periodic end

function create_cell_linesegments end
function create_cell_linesegments_periodic end

end
