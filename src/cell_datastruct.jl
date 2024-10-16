"Struct that stores indices of elements surrounding Cell (vertices, neighboring cells, and edges)"
struct CellConnectivity{MAX_N_EDGES,TI<:Integer}
    vertices::ImVecArray{MAX_N_EDGES, TI, 1}
    cells::ImVecArray{MAX_N_EDGES, TI, 1}
    edges::ImVecArray{MAX_N_EDGES, TI, 1}
end

function CellConnectivity(vertices::ImVecArray{N, TI, 1}) where  {N, TI}
    cells = ImmutableVectorArray(similar(vertices.data), vertices.length)
    edges = ImmutableVectorArray(similar(vertices.data), vertices.length)
    fill!(reinterpret(reshape,TI,cells.data), zero(TI))
    fill!(reinterpret(reshape,TI,edges.data), zero(TI))
    return CellConnectivity(vertices, cells, edges)
end

max_n_edges(::Type{<:CellConnectivity{N}}) where N = N
integer_precision(::Type{<:CellConnectivity{N,T}}) where {N,T} = T

Base.getproperty(cell::CellConnectivity,s::Symbol) = _getproperty(cell,Val(s))
_getproperty(cell::CellConnectivity,::Val{s}) where s = getfield(cell,s)

_getproperty(cell::CellConnectivity,::Val{:verticesOnCell}) = getfield(cell,:vertices)
_getproperty(cell::CellConnectivity,::Val{:edgesOnCell}) = getfield(cell,:edges)
_getproperty(cell::CellConnectivity,::Val{:cellsOnCell}) = getfield(cell,:cells)

const names_in_cell_connectivity = (:verticesOnCell, :edgesOnCell, :cellsOnCell)

"""
Struct with Cell's base information

# Fields

- `n::Integer`: Number of cells

- `indices::CellConnectivity`: Struct containing connectivitiy info

- `nEdges::Vector{UInt8}`: Vector with the number of edges (and vertices) of each cell

- `position::VecArray`: Vector with Cell's position points

- `periods::NTuple{2,Real}`: For periodic planar meshes, the (x,y) period, (0,0) if on the spehre

- `sphere_radius`: For spherical meshes, set to 0.0 if on the plane

- `onSphere`: Val{Bool} that tells if cells are on the sphere or not

"""
struct CellBase{S,MAX_N_EDGES,TI<:Integer,TF<:Real,Tz<:Number}
    n::Int
    indices::CellConnectivity{MAX_N_EDGES,TI} 
    nEdges::Vector{UInt8} 
    position::TensorsLite.VecMaybe2DxyArray{TF,Tz,1}
    periods::NTuple{2,TF}
    sphere_radius::TF
    onSphere::Val{S}

    function CellBase(position::Vec2DxyArray{TF,1}, indices::CellConnectivity{MN, TI}, xp::Number, yp::Number) where {TF, MN, TI}
        return new{false, MN, TI, TF, TensorsLite.Zeros.Zero}(length(position), indices, indices.vertices.length, position, (TF(xp) ,TF(yp)), TF(0.0), Val{false}())
    end

    function CellBase(position::Vec3DArray{TF,1}, indices::CellConnectivity{MN, TI}, r::Number) where {TF, MN, TI}
        return new{true, MN, TI, TF, TF}(length(position), indices, indices.vertices.length, position, (zero(TF), zero(TF)), r, Val{true}())
    end
end

on_a_sphere(::Type{<:CellBase{B}}) where B = B
max_n_edges(::Type{<:CellBase{B,N}}) where {N,B} = N
integer_precision(::Type{<:CellBase{N,B,T}}) where {N,B,T} = T
float_precision(::Type{<:CellBase{N,B,T,TF}}) where {N,B,T,TF} = TF

Base.getproperty(cell::CellBase,s::Symbol) = _getproperty(cell,Val(s))
_getproperty(cell::CellBase,::Val{s}) where s = getfield(cell,s)

for s in names_in_cell_connectivity
    @eval _getproperty(cell::CellBase,::Val{$(QuoteNode(s))}) = getproperty(getfield(cell,:indices),$(QuoteNode(s)))
end

_getproperty(cell::CellBase,::Val{:nCells}) = getfield(cell,:n)
_getproperty(cell::CellBase,::Val{:nEdgesOnCell}) = getfield(cell,:nEdges)
_getproperty(cell::CellBase,::Val{:xCell}) = getfield(cell,:position).x
_getproperty(cell::CellBase,::Val{:yCell}) = getfield(cell,:position).y
_getproperty(cell::CellBase,::Val{:zCell}) = getfield(cell,:position).z
_getproperty(cell::CellBase,::Val{:x_period}) = getfield(cell,:periods)[1]
_getproperty(cell::CellBase,::Val{:y_period}) = getfield(cell,:periods)[2]

const names_in_cell_base = (names_in_cell_connectivity..., :nCells, :nEdgesOnCell, :xCell, :yCell, :zCell, :x_period, :y_period)

"Struct holding Cell's information"
struct Cells{S,MAX_N_EDGES,TI<:Integer,TF<:Real,Tz<:Number}
    "Struct containing Cell's essential information"
    base::CellBase{S,MAX_N_EDGES,TI,TF,Tz}
    "Longitiude of Cell center (if on the sphere)"
    longitude::Vector{TF}
    "Latitude of Cell center (if on the sphere)"
    latitude::Vector{TF}
    """Cell's area"""
    area::Vector{TF}
    """Cartesian components of the vector pointing in the local vertical direction for a cell"""
    verticalUnitVectors::TensorsLite.VecMaybe2DxyArray{Tz,TF,1} # Actually, maybe 1Dz
    """Tuple with pair of Vectors of Vec3D's structs defining the tangent plane at a cell"""
    tangentPlane::NTuple{2,TensorsLite.VecMaybe2DxyArray{TF,Tz,1}}
    """Mesh density function evaluated at cell (used when generating the cell)"""
    meshDensity::Vector{TF}
end

on_a_sphere(::Type{<:Cells{B}}) where B = B
max_n_edges(::Type{<:Cells{B,N}}) where {B,N} = N
integer_precision(::Type{<:Cells{N,B,T}}) where {N,B,T} = T
float_precision(::Type{<:Cells{N,B,T,TF}}) where {N,B,T,TF} = TF

for fun in (:on_a_sphere, :max_n_edges, :integer_precision, :float_precision)
    @eval $fun(::T) where {T <: Union{<:Cells, <:CellBase}} = $fun(T)
end

Base.getproperty(cell::Cells,s::Symbol) = _getproperty(cell,Val(s))
_getproperty(cell::Cells,::Val{s}) where s = getfield(cell,s)

for s in (fieldnames(CellBase)..., names_in_cell_base...)
    @eval _getproperty(cell::Cells,::Val{$(QuoteNode(s))}) = _getproperty(getfield(cell,:base),Val{$(QuoteNode(s))}())
end

for (s,nc) in pairs((longitude=:lonCell, latitude=:latCell, area=:areaCell,
                     verticalUnitVectors=:localVerticalUnitVectors,
                     tangentPlane=:cellTangentPlane))

    @eval _getproperty(cell::Cells,::Val{$(QuoteNode(nc))}) = getfield(cell,$(QuoteNode(s)))
end

const names_in_cells = (names_in_cell_base..., :lonCell, :latCell, :meshDensity, :areaCell, :localVerticalUnitVectors, :cellTangentPlane)

#struct MPASCells{S,MAX_N_EDGES,TI<:Integer,TF<:Real,Tz<:Number}
#    base::Cells{S, MAX_N_EDGES, TI, TF, Tz}
#    """Mapping from local array index to global cell ID"""
#    indexToID::Vector{TI}
#    """Indicator of whether a cell is an interior cell, a relaxation-zone cell, or a specified-zone cell"""
#    bdyMask::Vector{TI}
#end
#
#Base.getproperty(cell::MPASCells,s::Symbol) = _getproperty(cell,Val(s))
#_getproperty(cell::MPASCells,::Val{s}) where s = getfield(cell,s)
#
#for s in (fieldnames(Cells)..., fieldnames(CellBase)...)
#    @eval _getproperty(cell::MPASCells,::Val{$(QuoteNode(s))}) = _getproperty(getfield(cell,:base),$(QuoteNode(s)))
#end
#
#_getproperty(cell::MPASCells, ::Val{:bdyMaskCell}) = getfield(cell, :bdyMask)
#_getproperty(cell::MPASCells, ::Val{:IndexToCellID}) = getfield(cell, :indexToID)
#
#const names_in_MPASCells = (names_in_cells..., :bdyMaskCell, :indexToCellID)
