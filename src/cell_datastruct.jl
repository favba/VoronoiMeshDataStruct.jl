struct CellConnectivity{MAX_N_EDGES,TI<:Integer}
    """Indices of vertices forming the cell"""
    vertices::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
    """Indices of the edges forming the cell"""
    edges::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
    """Indices of cells surrounding the given cell"""
    cells::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
end

max_n_edges(::Type{<:CellConnectivity{N}}) where N = N
integer_precision(::Type{<:CellConnectivity{N,T}}) where {N,T} = T

struct CellBase{MAX_N_EDGES, TI<:Integer, VAPos<:VecArray{<:Any,1}}
    """Cells connectivity data struct"""
    indices::CellConnectivity{MAX_N_EDGES,TI} 
    """Cell's number of edges (and vertices)"""
    nEdges::Vector{TI} 
    """Cell's x,y,z coordinates"""
    position::VAPos
end

max_n_edges(::Type{<:CellBase{N}}) where N = N
integer_precision(::Type{<:CellBase{N,T}}) where {N,T} = T
float_precision(::Type{<:CellBase{N,T,V}}) where {N,T,V} = TensorsLite._my_eltype(eltype(V))

mutable struct CellInfo{TCells<:CellBase,TI<:Integer,TF<:Real}
    const base::TCells
    longitude::Vector{TF}
    latitude::Vector{TF}
    """Mesh density function evaluated at cell (used when generating the cell)"""
    meshDensity::Vector{TF}
    """Mapping from local array index to global cell ID"""
    indexToID::Vector{TI}
    """Cell's area"""
    area::Vector{TF}
    """Indicator of whether a cell is an interior cell, a relaxation-zone cell, or a specified-zone cell"""
    bdyMask::Vector{TI}
    """Cartesian components of the vector pointing in the local vertical direction for a cell"""
    verticalUnitVectors::Vec3DArray{TF,1}
    """Tuple with pair of Vectors of Vec3D's structs defining the tangent plane at a cell"""
    tangentPlane::NTuple{2,Vec3DArray{TF,1}}
    """Coefficients for computing the off-diagonal components of the horizontal deformation"""
    defc_a::Matrix{TF}
    """Coefficients for computing the diagonal components of the horizontal deformation"""
    defc_b::Matrix{TF}
    """Coefficients for computing the x (zonal) derivative of a cell-centered variable"""
    xGradientCoeff::Matrix{TF}
    """Coefficients for computing the y (meridional) derivative of a cell-centered variable"""
    yGradientCoeff::Matrix{TF}
    """Coefficients to reconstruct velocity vectors at cell centers"""
    coeffsReconstruct::Vec3DArray{TF,2}

    function CellInfo(cell::TCells) where TCells<:CellBase
        return new{TCells,integer_precision(TCells),float_precision(TCells)}(cell)
    end
end

@inline Base.getproperty(cell::CellInfo,s::Symbol) = _getproperty(cell,Val(s))

@inline _getproperty(cell::CellInfo,::Val{s}) where s = getfield(cell,s)
@inline _getproperty(cell::CellInfo,::Val{:indices}) = getfield(cell,:base).indices
@inline _getproperty(cell::CellInfo,::Val{:nEdges}) = getfield(cell,:base).nEdges
@inline _getproperty(cell::CellInfo,::Val{:position}) = getfield(cell,:base).position