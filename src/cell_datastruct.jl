struct CellConnectivity{MAX_N_EDGES,TI<:Integer}
    """Indices of vertices forming the cell"""
    vertices::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
    """Indices of the edges forming the cell"""
    edges::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
    """Indices of cells surrounding the given cell"""
    cells::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
end

Base.getproperty(cell::CellConnectivity,s::Symbol) = _getproperty(cell,Val(s))
_getproperty(cell::CellConnectivity,::Val{s}) where s = getfield(cell,s)
_getproperty(cell::CellConnectivity,::Val{:verticesOnCell}) = getfield(cell,:vertices)
_getproperty(cell::CellConnectivity,::Val{:edgesOnCell}) = getfield(cell,:edges)
_getproperty(cell::CellConnectivity,::Val{:cellsOnCell}) = getfield(cell,:cells)

max_n_edges(::Type{<:CellConnectivity{N}}) where N = N
integer_precision(::Type{<:CellConnectivity{N,T}}) where {N,T} = T

struct CellBase{MAX_N_EDGES, TI<:Integer, VAPos<:VecArray{<:Any,1}}
    n::Int
    """Cells connectivity data struct"""
    indices::CellConnectivity{MAX_N_EDGES,TI} 
    """Cell's number of edges (and vertices)"""
    nEdges::Vector{TI} 
    """Cell's x,y,z coordinates"""
    position::VAPos
end

Base.getproperty(cell::CellBase,s::Symbol) = _getproperty(cell,Val(s))
_getproperty(cell::CellBase,::Val{s}) where s = getfield(cell,s)
for s in (:verticesOnCell,:edgesOnCell,:cellsOnCell)
    @eval _getproperty(cell::CellBase,::Val{$(QuoteNode(s))}) = getproperty(getfield(cell,:indices),$(QuoteNode(s)))
end
_getproperty(cell::CellBase,::Val{:nCells}) = getfield(cell,:n)
_getproperty(cell::CellBase,::Val{:nEdgesOnCell}) = getfield(cell,:nEdges)
_getproperty(cell::CellBase,::Val{:xCell}) = getfield(cell,:position).x
_getproperty(cell::CellBase,::Val{:yCell}) = getfield(cell,:position).y
_getproperty(cell::CellBase,::Val{:zCell}) = getfield(cell,:position).z

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
    defcA::Matrix{TF}
    """Coefficients for computing the diagonal components of the horizontal deformation"""
    defcB::Matrix{TF}
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

Base.getproperty(cell::CellInfo,s::Symbol) = _getproperty(cell,Val(s))

_getproperty(cell::CellInfo,::Val{s}) where s = getfield(cell,s)
_getproperty(cell::CellInfo,::Val{:n}) = getfield(cell,:base).n
_getproperty(cell::CellInfo,::Val{:indices}) = getfield(cell,:base).indices
_getproperty(cell::CellInfo,::Val{:nEdges}) = getfield(cell,:base).nEdges
_getproperty(cell::CellInfo,::Val{:position}) = getfield(cell,:base).position

for s in (:nCells,:verticesOnCell,:edgesOnCell,:cellsOnCell,:xCell,:yCell,:zCell,:nEdgesOnCell)
    @eval _getproperty(cell::CellInfo,::Val{$(QuoteNode(s))}) = getproperty(getfield(cell,:base),$(QuoteNode(s)))
end

for (s,nc) in pairs((longitude=:lonCell, latitude=:latCell, meshDensity=:meshDensity,
                     indexToID=:indexToCellID, area=:areaCell, bdyMask=:bdyMaskCell,
                     verticalUnitVectors=:localVerticalUnitVectors, tangentPlane=:cellTangentPlane,
                     defcA=:defc_a, defcB=:defc_b, xGradientCoeff=:cell_gradient_coef_x,
                     yGradientCoeff=:cell_gradient_coef_y, coeffsReconstruct=:coeffs_reconstruct))

    @eval _getproperty(cell::CellInfo,::Val{$(QuoteNode(nc))}) = getfield(cell,$(QuoteNode(s)))
end