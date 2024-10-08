"Struct that stores indices of elements surrounding Cell (vertices, edges and neighboring cells)"
struct CellConnectivity{MAX_N_EDGES,TI<:Integer}
    """Indices of vertices forming the cell"""
    vertices::Vector{ImmutableVector{MAX_N_EDGES,TI}}
    """Indices of the edges forming the cell"""
    edges::Vector{ImmutableVector{MAX_N_EDGES,TI}}
    """Indices of cells surrounding the given cell"""
    cells::Vector{ImmutableVector{MAX_N_EDGES,TI}}
end

Base.getproperty(cell::CellConnectivity,s::Symbol) = _getproperty(cell,Val(s))
_getproperty(cell::CellConnectivity,::Val{s}) where s = getfield(cell,s)
_getproperty(cell::CellConnectivity,::Val{:verticesOnCell}) = getfield(cell,:vertices)
_getproperty(cell::CellConnectivity,::Val{:edgesOnCell}) = getfield(cell,:edges)
_getproperty(cell::CellConnectivity,::Val{:cellsOnCell}) = getfield(cell,:cells)

max_n_edges(::Type{<:CellConnectivity{N}}) where N = N
integer_precision(::Type{<:CellConnectivity{N,T}}) where {N,T} = T

"Struct with Cell's essential information"
struct CellBase{S,MAX_N_EDGES,TI<:Integer,TF<:Real,Tz<:Number}
    "Number of Cells"
    n::Int
    """Cells connectivity data struct"""
    indices::CellConnectivity{MAX_N_EDGES,TI} 
    """Cell's number of edges (and vertices)"""
    nEdges::Vector{TI} 
    """Cell's x,y,z coordinates"""
    position::TensorsLite.VecMaybe2DxyArray{TF,Tz,1}
    "Val{Bool} that tells if cells are on the sphere or not"
    onSphere::Val{S}
end

on_a_sphere(::Type{<:CellBase{B}}) where B = B
max_n_edges(::Type{<:CellBase{B,N}}) where {N,B} = N
integer_precision(::Type{<:CellBase{N,B,T}}) where {N,B,T} = T
float_precision(::Type{<:CellBase{N,B,T,TF}}) where {N,B,T,TF} = TF

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

"Struct holding all Cell's information (may be incompletely initialized)"
mutable struct CellInfo{S,MAX_N_EDGES,TI<:Integer,TF<:Real,Tz<:Number}
    "Struct containing Cell's essential information"
    const base::CellBase{S,MAX_N_EDGES,TI,TF,Tz}
    "Longitiude of Cell center (if on the sphere)"
    longitude::Vector{TF}
    "Latitude of Cell center (if on the sphere)"
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
    verticalUnitVectors::TensorsLite.VecMaybe2DxyArray{Tz,TF,1} # Actually, maybe 1Dz
    """Tuple with pair of Vectors of Vec3D's structs defining the tangent plane at a cell"""
    tangentPlane::NTuple{2,TensorsLite.VecMaybe2DxyArray{TF,Tz,1}}
    """Coefficients for computing the off-diagonal components of the horizontal deformation"""
    defcA::Matrix{TF}
    """Coefficients for computing the diagonal components of the horizontal deformation"""
    defcB::Matrix{TF}
    """Coefficients for computing the x (zonal) derivative of a cell-centered variable"""
    xGradientCoeff::Matrix{TF}
    """Coefficients for computing the y (meridional) derivative of a cell-centered variable"""
    yGradientCoeff::Matrix{TF}
    """Coefficients to reconstruct velocity vectors at cell centers"""
    coeffsReconstruct::TensorsLite.VecMaybe2DxyArray{TF,Tz,2}

    function CellInfo(cell::CellBase{N,S,TI,TF,Tz}) where {N,S,TI,TF,Tz}
        return new{N,S,TI,TF,Tz}(cell)
    end
end

on_a_sphere(::Type{<:CellInfo{B}}) where B = B
max_n_edges(::Type{<:CellInfo{B,N}}) where {B,N} = N
integer_precision(::Type{<:CellInfo{N,B,T}}) where {N,B,T} = T
float_precision(::Type{<:CellInfo{N,B,T,TF}}) where {N,B,T,TF} = TF

for fun in (:on_a_sphere, :max_n_edges, :integer_precision, :float_precision)
    @eval $fun(::T) where {T <: Union{<:CellInfo, <:CellBase}} = $fun(T)
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
