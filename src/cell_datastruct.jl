struct CellConnectivity{MAX_N_EDGES,TI<:Integer}
    """Indices of vertices forming the cell"""
    vertices::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
    """Indices of the edges forming the cell"""
    edges::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
    """Indices of cells surrounding the given cell"""
    cells::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
end

max_n_edges(::CellConnectivity{N}) where N = N
integer_precision(::CellConnectivity{N,T}) where {N,T} = T

struct Cells{MAX_N_EDGES, TI<:Integer, VAPos<:VecArray{<:Any,1}}
    """Cells connectivity data struct"""
    indices::CellConnectivity{MAX_N_EDGES,TI} 
    """Cell's number of edges (and vertices)"""
    nEdges::Vector{TI} 
    """Cell's x,y,z coordinates"""
    position::VAPos
end

#struct CellDataStruct{MAX_N_EDGES, TF<:AbstractFloat, TI<:Integer}
#    """Cell's x,y,z"""
#    xyzCell::Vec3DArray{TF,1}
#    """Cell's longitude and latitude coordinates encapsulated in a Vec2Dxy struct"""
#    lonlatCell::Vec2DxyArray{TF,1}
#    """Mapping from local array index to global cell ID"""
#    indexToCellID::Vector{TI} 
#
#    """Cell's number of edges"""
#    nEdgesOnCell::Vector{TI} 
#    """Indices of the edges forming the cell"""
#    edgesOnCell::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
#    """Cell's area"""
#    areaCell::Vector{TF}
#    """Indices of cells surrounding the given cell"""
#    cellsOnCell::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
#    """Indices of vertices forming the cell"""
#    verticesOnCell::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
#    """Mesh density function evaluated at cell (used when generating the cell)"""
#    meshDensity::Vector{TF}
#    #"""Indicator of whether a cell is an interior cell, a relaxation-zone cell, or a specified-zone cell"""
#    #bdyMaskCell::Vector{TI}
#    """Cartesian components of the vector pointing in the local vertical direction for a cell"""
#    localVerticalUnitVectors::Vec3DArray{TF,1}
#    """Tuple with pair of Vectors of Vec3D's structs defining the tangent plane at a cell"""
#    cellTangentPlane::NTuple{2,Vec3DArray{TF,1}}
#end